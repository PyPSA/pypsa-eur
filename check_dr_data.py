"""
check_dr_data.py
================
Checks that prepare_dr_data.py found, read and wrote everything the config
asks for. Run it after prepare_dr_data.py, with the same config:

    python check_dr_data.py                 # uses config_dr.yaml
    python check_dr_data.py --config my.yaml

It checks:
  1. Raw files    - every needed Zenodo file is in data_dir and has the
                    expected layout (regions x 12 months x 24 hours, years).
  2. Output files - every expected output file exists.
  3. Contents     - right number of hours, no missing or negative values,
                    the right regions, every day in a month identical.
  4. Values       - random hours compared with the raw Zenodo numbers.

Ends with PASS or FAIL. Needs no internet and no PyPSA.
"""

import argparse
import os
import random
import re
import sys

import pandas as pd
import yaml

FILE_PATTERN = re.compile(r"^([dp])_([A-Z]{2})V(\d+)\.csv$")
ISO2_TO_EUROSTAT = {"GR": "EL", "GB": "UK"}
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

problems = []


def ok(msg):
    print(f"  ok    {msg}")


def fail(msg):
    print(f"  FAIL  {msg}")
    problems.append(msg)


def find_file(data_dir, prefix, device):
    hits = []
    for name in os.listdir(data_dir):
        m = FILE_PATTERN.match(name)
        if m and m.group(1) == prefix and m.group(2) == device:
            hits.append((int(m.group(3)), name))
    return os.path.join(data_dir, max(hits)[1]) if hits else None



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config_dr.yaml")
    ap.add_argument("--samples", type=int, default=20,
                    help="random values to compare per output file")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    countries = ["NO" if c is False else str(c).upper()
                 for c in (cfg.get("countries") or [])]
    codes = {ISO2_TO_EUROSTAT.get(c, c) for c in countries}
    data_dir = cfg["retrieve"]["data_dir"]
    variables = cfg["variables"]
    need_p = bool({"capacity", "max_increase"} & set(variables))
    use_rates = bool((cfg.get("participation") or {}).get("apply"))
    out_dir = os.path.join(cfg["output"]["dir"],
                           "achievable" if use_rates else "theoretical")
                           
    n_hours = 366 * 24  # month, day, hour incl. 29 February

    print(f"Config: {args.config}")
    print(f"  devices={cfg['devices']} years={cfg['years']} variables={variables}")
    print(f"  countries={countries or 'all'} ({n_hours} rows: month, day, hour)")
    print(f"  raw data: {data_dir}   output: {out_dir}\n")

    # participation rates (only if used)
    rates = None
    if use_rates:
        pcfg = cfg["participation"]
        pr = pd.read_csv(os.path.join(data_dir, pcfg["file"]), index_col=0)
        pr = pr.drop(index="Average", errors="ignore")
        pr.index = pr.index.map(COUNTRY_NAMES)
        rates = {d: pr[pcfg["device_groups"][d]] / 100 for d in cfg["devices"]}

    # ---- 1. raw files -----------------------------------------------------
    print("1. Raw Zenodo files")
    raw = {}
    for extra in cfg["retrieve"].get("extra_files") or []:
        if os.path.exists(os.path.join(data_dir, extra)):
            ok(extra)
        else:
            fail(f"{extra} missing in {data_dir}")

    for dev in cfg["devices"]:
        for prefix in ["d"] + (["p"] if need_p else []):
            path = find_file(data_dir, prefix, dev)
            if path is None:
                fail(f"{prefix}_{dev}V*.csv missing in {data_dir}")
                continue
            df = pd.read_csv(path, index_col=0)
            name = os.path.basename(path)
            years_missing = [y for y in cfg["years"] if str(y) not in df.columns]
            n_reg = df["nutscode"].nunique()
            bad_layout = (df.groupby(["nutscode", "TIME"]).size() != 24).any() \
                or df["TIME"].nunique() != 12
            if years_missing:
                fail(f"{name}: years {years_missing} not in file")
            elif bad_layout:
                fail(f"{name}: not 12 months x 24 hours per region")
            else:
                ok(f"{name}: {n_reg} regions x 12 months x 24 hours, years present")
            df["month"] = df["TIME"].str[-2:].astype(int)
            raw[(prefix, dev)] = df

    # ---- 2-4. outputs -----------------------------------------------------
    print("\n2-4. Output files")
    for dev in cfg["devices"]:
        if ("d", dev) not in raw:
            continue
        d = raw[("d", dev)]
        if codes:
            d = d[d["country"].isin(codes)]
        regions = sorted(d["nutscode"].unique())

        for var in variables:
            if var != "load" and ("p", dev) not in raw:
                continue
            for y in cfg["years"]:
                fn = os.path.join(out_dir, f"{var}_{dev}_{y}.csv")
                if not os.path.exists(fn):
                    fail(f"{os.path.basename(fn)} missing")
                    continue
                f = pd.read_csv(fn, index_col=["month", "day", "hour"])
                name = os.path.basename(fn)
                errs = []

                # 3. contents
                if len(f) != n_hours:
                    errs.append(f"{len(f)} rows, expected {n_hours}")
                if f.isna().any().any():
                    errs.append(f"{int(f.isna().sum().sum())} missing values")
                if (f < 0).any().any():
                    errs.append("negative values")
                if cfg.get("nuts_to_bus") is None and sorted(f.columns) != regions:
                    errs.append(f"{f.shape[1]} columns, expected {len(regions)} regions")
                for m, grp in f.groupby(level="month"):
                    days = grp.to_numpy().reshape(-1, 24, grp.shape[1])
                    if not (days == days[0]).all():
                        errs.append(f"month {m}: days are not identical")
                        break

                # 4. values vs raw data (only for NUTS2 output)
                if cfg.get("nuts_to_bus") is None and not errs:
                    src = {"load": d}
                    if var != "load":
                        p = raw[("p", dev)]
                        src["cap"] = p[p["country"].isin(codes)] if codes else p
                    dec = cfg["output"].get("decimals", 4)
                    random.seed(0)
                    for _ in range(args.samples):
                        ts = random.choice(f.index)
                        reg = random.choice(regions)

                        def raw_val(df):
                            row = df[(df.nutscode == reg) & (df.month == ts[0]) & (df.hour == ts[2])]
                            v = row[str(y)].iloc[0]
                            if pd.isna(v):  # the script fills gaps from nearby years
                                ycols = [c for c in df.columns if c.isdigit()]
                                v = row[ycols].ffill(axis=1).bfill(axis=1)[str(y)].iloc[0]
                            return v

                        if var == "load":
                            expect = raw_val(src["load"])
                        else:
                            cap = raw_val(src["cap"])
                            expect = cap if var == "capacity" else max(cap - raw_val(src["load"]), 0)
                        if rates is not None:
                            expect *= rates[dev][reg[:2]]
                        got = f.loc[ts, reg]
                        if abs(got - expect) > 10 ** -dec:
                            errs.append(f"{reg} {ts}: file {got}, raw data {expect:.{dec}f}")
                            break

                if errs:
                    fail(f"{name}: " + "; ".join(errs))
                else:
                    ok(f"{name}: {len(f)} hours x {f.shape[1]} regions, "
                       f"{args.samples} values match raw data")

    print("\n" + ("PASS - everything found, read and written correctly"
                  if not problems else f"FAIL - {len(problems)} problem(s) listed above"))
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
