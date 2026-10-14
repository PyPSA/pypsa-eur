"""
Retrieve and prepare the O'Reilly et al. (2022) openENTRANCE residential
demand-response dataset (Zenodo 10.5281/zenodo.7871106) for PyPSA-Eur.
Everything is controlled by a YAML config file (default: config_dr.yaml).

Steps
-----
1. Retrieve: download the d_/p_ files for the configured devices from Zenodo
   (only files that are missing or fail the checksum), using the Zenodo API.
2. Select:   keep the configured years and countries.
3. Expand:   each month has one representative 24-hour day; repeat it for
             every day of that month to get a full hourly year that matches
             the PyPSA-Eur snapshots.
4. Scale (optional): multiply by country participation rates to go from
             theoretical to achievable potential.
5. Aggregate (optional): NUTS2 regions -> PyPSA-Eur buses.
6. Write:    one CSV per variable, device and year:
                 <output.dir>/{theoretical|achievable}/<variable>_<DEVICE>_<year>.csv
             rows = hourly snapshots, columns = NUTS2 regions (or buses)

"""

import argparse
import hashlib
import os
import re
import sys
import time

import pandas as pd
import requests
import yaml

DEVICES = {
    "AC": "air conditioning",
    "CP": "heat circulation pump",
    "DW": "dish washer",
    "EV": "battery electric vehicle",
    "FR": "freezer",
    "HP": "air-to-air heat pump",
    "RF": "refrigerator",
    "SH": "storage heater",
    "TD": "tumble drier",
    "WH": "water heater",
    "WM": "washing machine",
}

# The dataset uses Eurostat codes and PyPSA-Eur uses ISO2, this translates between them: 
EUROSTAT_TO_ISO2 = {"EL": "GR", "UK": "GB"}
ISO2_TO_EUROSTAT = {v: k for k, v in EUROSTAT_TO_ISO2.items()}

# Country names in participation_rates_country.csv translated to Eurostat codes
COUNTRY_NAMES = {
    "Austria": "AT", "Belgium": "BE", "Bulgaria": "BG", "Croatia": "HR",
    "Cyprus": "CY", "Czechia": "CZ", "Denmark": "DK", "Estonia": "EE",
    "Finland": "FI", "France": "FR", "Germany": "DE", "Greece": "EL",
    "Hungary": "HU", "Ireland": "IE", "Italy": "IT", "Latvia": "LV",
    "Lithuania": "LT", "Luxembourg": "LU", "Malta": "MT", "Norway": "NO",
    "Poland": "PL", "Portugal": "PT", "Romania": "RO", "Slovakia": "SK",
    "Slovenia": "SI", "Spain": "ES", "Sweden": "SE", "Switzerland": "CH",
    "Netherlands": "NL", "Turkey": "TR", "United Kingdom": "UK",
}

FILE_PATTERN = re.compile(r"^([dp])_([A-Z]{2})V(\d+)\.csv$")



# 1. Retrieve

def zenodo_file_list(record_id):
    """Return {filename: (download_url, md5)} for a Zenodo record."""
    url = f"https://zenodo.org/api/records/{record_id}"
    for attempt in range(5):
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            break
        except requests.RequestException as e:
            wait = 10 * (attempt + 1)
            print(f"  Zenodo API not reachable ({e}); retrying in {wait}s")
            time.sleep(wait)
    else:
        raise RuntimeError(f"Could not reach {url}")

    files = {}
    for f in r.json()["files"]:
        md5 = f.get("checksum", "").removeprefix("md5:")
        link = f["links"].get("self") or f["links"].get("content")
        files[f["key"]] = (link, md5)
    return files


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url, dest, md5=None, retries=5):
    tmp = dest + ".part"
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                done = 0
                with open(tmp, "wb") as fh:
                    for block in r.iter_content(1 << 20):
                        fh.write(block)
                        done += len(block)
                        if total:
                            print(f"\r    {done / 1e6:6.1f} / {total / 1e6:.1f} MB", end="")
                print()
            if md5 and md5sum(tmp) != md5:
                raise IOError("checksum mismatch")
            os.replace(tmp, dest)
            return
        except (requests.RequestException, IOError) as e:
            wait = 15 * (attempt + 1)
            print(f"    failed ({e}); retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"Could not download {url}")


def retrieve(cfg):
    rcfg = cfg["retrieve"]
    data_dir = rcfg["data_dir"]
    os.makedirs(data_dir, exist_ok=True)

    print(f"Retrieving file list for Zenodo record {rcfg['zenodo_record']}")
    available = zenodo_file_list(rcfg["zenodo_record"])

    wanted = set(rcfg.get("extra_files") or [])
    prefixes = {"d"} | ({"p"} if needs_capacity(cfg) else set())
    for name in available:
        m = FILE_PATTERN.match(name)
        if m and m.group(1) in prefixes and m.group(2) in cfg["devices"]:
            wanted.add(name)

    for name in sorted(wanted):
        if name not in available:
            print(f"  ! {name} is not in the Zenodo record")
            continue
        url, md5 = available[name]
        dest = os.path.join(data_dir, name)
        if os.path.exists(dest) and (not md5 or md5sum(dest) == md5):
            print(f"  ok   {name}")
            continue
        print(f"  get  {name}")
        download(url, dest, md5)


# 2-4. Process

def needs_capacity(cfg):
    return bool({"capacity", "max_increase"} & set(cfg["variables"]))


def find_file(data_dir, prefix, device):
    """Find e.g. d_ACV8.csv regardless of version; take the highest version."""
    hits = []
    for name in os.listdir(data_dir):
        m = FILE_PATTERN.match(name)
        if m and m.group(1) == prefix and m.group(2) == device:
            hits.append((int(m.group(3)), name))
    return os.path.join(data_dir, max(hits)[1]) if hits else None


def read_representative_days(path, years, countries):
    """
    Return {year: DataFrame} with index (month, hour) = 288 rows and one
    column per NUTS2 region.
    """
    df = pd.read_csv(path, index_col=0)

    missing = [y for y in years if str(y) not in df.columns]
    if missing:
        raise ValueError(f"{os.path.basename(path)}: years {missing} not in file")

    if countries:
        codes = {ISO2_TO_EUROSTAT.get(c, c) for c in countries}
        df = df[df["country"].isin(codes)]
        if df.empty:
            raise ValueError(f"No rows left after filtering countries {countries}")

    df["month"] = df["TIME"].str[-2:].astype(int)  # '2018M01' -> 1
    df["hour"] = df["hour"].astype(int)

    # Gaps along the year axis (e.g. CP / EL30 has no 2049-2050):
    # carry the nearest available year forward, then backward.
    year_cols = [c for c in df.columns if c.isdigit()]
    n_nan = int(df[year_cols].isna().sum().sum())
    df[year_cols] = df[year_cols].ffill(axis=1).bfill(axis=1)
    if n_nan:
        print(f"    filled {n_nan} missing values in {os.path.basename(path)}")

    out = {}
    for y in years:
        wide = df.pivot_table(index=["month", "hour"], columns="nutscode",
                              values=str(y), aggfunc="first")
        if wide.shape[0] != 12 * 24:
            raise ValueError(f"{path}: expected 288 (month, hour) rows, got {wide.shape[0]}")
        out[y] = wide.fillna(0.0)
    return out


def expand_to_year(rep_days):
    """
    Repeat each month's representative day for every day of that month.
    Rows are labelled (month, day, hour), no year. 29 February is included
    (same as any other February day), so the data fits any year.
    """
    days = pd.date_range("2024-01-01", "2024-12-31 23:00", freq="h")  # 2024 = leap year
    # Any leap year works, it is only used to list all days (incl. 29 Feb). 
    full = rep_days.reindex(pd.MultiIndex.from_arrays([days.month, days.hour]))
    full.index = pd.MultiIndex.from_arrays([days.month, days.day, days.hour],
                                           names=["month", "day", "hour"])
    return full


def aggregate_to_buses(df, mapping):
    """mapping: columns nutscode, bus, optional weight (share of region per bus)."""
    m = mapping.copy()
    if "weight" not in m.columns:
        m["weight"] = 1.0
    m = m[m["nutscode"].isin(df.columns)]
    unmapped = sorted(set(df.columns) - set(m["nutscode"]))
    if unmapped:
        print(f"    {len(unmapped)} NUTS2 regions have no bus and are dropped: "
              f"{unmapped[:8]}{'...' if len(unmapped) > 8 else ''}")
    weights = m.pivot_table(index="nutscode", columns="bus", values="weight",
                            aggfunc="sum", fill_value=0.0)
    return df[weights.index] @ weights


def read_participation_rates(cfg):
    """
    Return {device: Series(rate as fraction, index = NUTS2 country code)}.
    The file has one row per country and one column per device group (%).
    """
    pcfg = cfg["participation"]
    path = os.path.join(cfg["retrieve"]["data_dir"], pcfg["file"])
    pr = pd.read_csv(path, index_col=0)
    pr = pr.drop(index="Average", errors="ignore")

    unknown = sorted(set(pr.index) - set(COUNTRY_NAMES))
    if unknown:
        raise ValueError(f"Unknown country names in {pcfg['file']}: {unknown}")
    pr.index = pr.index.map(COUNTRY_NAMES)

    rates = {}
    for dev in cfg["devices"]:
        group = pcfg["device_groups"].get(dev)
        if group not in pr.columns:
            raise ValueError(f"No participation column '{group}' for device {dev}")
        rates[dev] = pr[group] / 100.0
    return rates


def apply_rates(rep, rate):
    """Scale each NUTS2 column by its country's participation rate."""
    country = rep.columns.str[:2]
    missing = sorted(set(country) - set(rate.index))
    if missing:
        raise ValueError(f"No participation rate for countries {missing}")
    return rep * rate.reindex(country).to_numpy()


def process(cfg):
    data_dir = cfg["retrieve"]["data_dir"]
    pcfg = cfg.get("participation") or {}
    use_rates = bool(pcfg.get("apply"))
    rates = read_participation_rates(cfg) if use_rates else None
    out_dir = os.path.join(cfg["output"]["dir"],
                           "achievable" if use_rates else "theoretical")
    decimals = cfg["output"].get("decimals", 4)
    years = cfg["years"]
    countries = cfg.get("countries") or []
    variables = cfg["variables"]
    mapping = pd.read_csv(cfg["nuts_to_bus"]) if cfg.get("nuts_to_bus") else None
    os.makedirs(out_dir, exist_ok=True)

    written, skipped = [], []
    for dev in cfg["devices"]:
        d_path = find_file(data_dir, "d", dev)
        if d_path is None:
            skipped.append(f"{dev} (no d_ file)")
            continue
        p_path = find_file(data_dir, "p", dev) if needs_capacity(cfg) else None
        print(f"[{dev}] {DEVICES.get(dev, dev)}")

        load = read_representative_days(d_path, years, countries)
        cap = read_representative_days(p_path, years, countries) if p_path else None
        if needs_capacity(cfg) and cap is None:
            skipped.append(f"{dev} capacity/max_increase (no p_ file)")

        for y in years:
            series = {}
            if "load" in variables:
                series["load"] = load[y]
            if cap is not None:
                c = cap[y].reindex(columns=load[y].columns).fillna(0.0)
                if "capacity" in variables:
                    series["capacity"] = c
                if "max_increase" in variables:
                    series["max_increase"] = (c - load[y]).clip(lower=0.0)

            for var, rep in series.items():
                if use_rates:
                    rep = apply_rates(rep, rates[dev])
                full = expand_to_year(rep)
                if mapping is not None:
                    full = aggregate_to_buses(full, mapping)
                fn = os.path.join(out_dir, f"{var}_{dev}_{y}.csv")
                full.round(decimals).to_csv(fn)
                written.append(fn)

    print(f"\nWrote {len(written)} files to {out_dir}")
    if skipped:
        print("Skipped: " + "; ".join(skipped))
    return written







def load_config(path):
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    unknown = set(cfg["devices"]) - set(DEVICES)
    if unknown:
        raise ValueError(f"Unknown device codes in config: {sorted(unknown)}")
    # YAML reads an unquoted NO (Norway) as boolean false - turn it back
    cfg["countries"] = ["NO" if c is False else str(c).upper()
                        for c in (cfg.get("countries") or [])]
    bad_vars = set(cfg["variables"]) - {"load", "capacity", "max_increase"}
    if bad_vars:
        raise ValueError(f"Unknown variables in config: {sorted(bad_vars)}")
    return cfg


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config_dr.yaml")
    ap.add_argument("--retrieve-only", action="store_true")
    ap.add_argument("--no-retrieve", action="store_true",
                    help="skip download even if enabled in the config")
    args = ap.parse_args()

    cfg = load_config(args.config)
    if cfg["retrieve"].get("enable", True) and not args.no_retrieve:
        retrieve(cfg)
    if args.retrieve_only:
        return 0
    return 0 if process(cfg) else 1


if __name__ == "__main__":
    sys.exit(main())
