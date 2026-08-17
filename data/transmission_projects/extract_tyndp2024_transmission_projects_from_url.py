"""
Extract TYNDP 2024 transmission projects from the ENTSO-E project map
and the TYNDP 2024 Trans.Investments Excel sheet.

The resulting CSV files follow the format expected by PyPSA-Eur:

    data/transmission_projects/<project>/new_lines.csv
    data/transmission_projects/<project>/new_links.csv

The map API is used for project coordinates.
The Excel file is used for project metadata, including route length.
"""

from pathlib import Path
import requests
import pandas as pd


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MAP_URL = "https://tyndp2024.entsoe.eu/api/maps/all"

OUTPUT_DIR = Path("data/transmission_projects/tyndp2024")

# Change this to the location of your TYNDP 2024 Excel file
EXCEL_FILE = Path("data/transmission_projects/20250312_export_transmission.xlsx")

EXCEL_SHEET = "Trans.Investments"


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

STATUS_COLUMN = (
    "Status ID\n"
    "1 : Under Consideration,\n"
    "2 : In Planning but not permitting,\n"
    "3 : In permitting,\n"
    "4 : Under Construction\n"
    "5 : Commissioned\n"
    "6 : Completed"
)

TRANSMISSION_TYPES = [
    "ACTransmissionLine",
    "DCTransmissionLine",
    "OffshoreACTransmissionCable",
    "OffshoreDCTransmissionCable",
    "OnshoreACTransmissionCable"
]

def download_map():
    """Download the TYNDP 2024 project map."""

    print(f"Downloading project map from:\n{MAP_URL}")

    response = requests.get(MAP_URL, timeout=60)
    response.raise_for_status()

    data = response.json()

    if data.get("type") != "FeatureCollection":
        raise ValueError("Unexpected map format.")

    print(f"Downloaded {len(data['features'])} map features.\n")

    return data

def read_investments(excel_file):
    """Read the TYNDP 2024 Trans.Investments sheet."""

    print(f"\nReading Excel file:\n{excel_file}")

    investments = pd.read_excel(
        excel_file,
        sheet_name=EXCEL_SHEET,
        header=1,
    )

    print(f"Read {len(investments)} investments.")

    return investments

def create_map_dataframe(map_data):
    """Convert GeoJSON features into a dataframe."""

    rows = []

    for feature in map_data["features"]:

        properties = feature.get("properties", {})
        geometry = feature.get("geometry", {})

        row = properties.copy()

        row["geometry_type"] = geometry.get("type")
        row["coordinates"] = geometry.get("coordinates")

        rows.append(row)

    map_df = pd.DataFrame(rows)

    print(f"\nMap dataframe: {map_df.shape}")

    if "type" in map_df.columns:
        print("Map element types:")
        print(map_df["type"].value_counts(dropna=False))

    return map_df

def print_test_result(name, mask, dataframe):
    """Print how many rows pass/fail a test."""

    failed = (~mask).sum()
    passed = mask.sum()

    print(
        f"\nThe outcome of the {name} test. "
        f"Passed: {passed:>4}   "
        f"Failed: {failed:>4}"
    )

    if failed > 0:
        failed_rows = dataframe.loc[~mask].copy()

        failed_rows["Investment Name"] = (
            failed_rows["Investment Name"]
            .astype(str)
            .str.slice(0, 40)
        )

        failed_rows["Type of Element"] = (
            failed_rows["Type of Element"]
            .astype(str)
            .str.slice(0, 40)
        )

        print(
            failed_rows[
                [
                    "This investment belongs to project number…",
                    "Investment number",
                    "Investment Name",
                    "Type of Element",
                ]
            ].to_string(index=False)
        )

def create_projects(investments, map_df, max_build_year=2030):
    """
    Select transmission investments and combine Excel metadata
    with map coordinates.
    """

    print("\n" + "=" * 80)
    print("CREATING TRANSMISSION PROJECTS")
    print("=" * 80)

    print(f"\nInput investments: {len(investments)}")

    # ------------------------------------------------------------------
    # 1. Transmission elements only
    # ------------------------------------------------------------------

    mask = investments["Type of Element"].isin(TRANSMISSION_TYPES)

    print_test_result(
        "transmission element type",
        mask,
        investments,
        )

    trans = investments.loc[mask].copy()

    # ------------------------------------------------------------------
    # 2. Investment ID must exist
    # ------------------------------------------------------------------

    mask = trans["Investment number"].notna()

    print_test_result(
        "investment number exists",
        mask,
        trans,
        )

    trans = trans.loc[mask].copy()

    # ------------------------------------------------------------------
    # 3. Map feature must exist
    # ------------------------------------------------------------------

    map_ids = set(
        map_df["investmentId"]
        .dropna()
        .astype(int)
        )

    mask = trans["Investment number"].astype(int).isin(map_ids)

    print_test_result(
        "investment exists on project map",
        mask,
        trans,
        )

    trans = trans.loc[mask].copy()

    # ------------------------------------------------------------------
    # 4. Geometry must exist
    # ------------------------------------------------------------------

    map_lookup = (
        map_df
        .dropna(subset=["investmentId"])
        .copy()
        )

    map_lookup["investmentId"] = (
        map_lookup["investmentId"]
        .astype(int)
        )

    map_lookup = map_lookup.drop_duplicates(
        subset="investmentId"
        )

    trans = trans.merge(
        map_lookup[
            [
                "investmentId",
                "geometry_type",
                "coordinates",
            ]
        ],
        left_on="Investment number",
        right_on="investmentId",
        how="left",
        )

    mask = trans["coordinates"].notna()

    print_test_result(
        "map geometry exists",
        mask,
        trans,
        )

    trans = trans.loc[mask].copy()

    # ------------------------------------------------------------------
    # 5. Commissioning year
    #
    # We introduce this check because the further away the projects
    # are, the more uncertain the actual project geometry is, and thus
    # it is not easy to create a line project for this.
    # ------------------------------------------------------------------

    commissioning_dates = pd.to_datetime(
        trans["Commissioning Year"],
        format="%m-%Y",
        errors="coerce",
        )

    mask = (
        commissioning_dates.notna()
        & (commissioning_dates.dt.year <= max_build_year)
        )

    print_test_result(
        f"commissioned by {max_build_year}",
        mask,
        trans,
        )

    trans = trans.loc[mask].copy()

    # ------------------------------------------------------------------
    # 6. Route length must exist
    # ------------------------------------------------------------------

    mask = trans["Total route length (km)"].notna()

    print_test_result(
        "route length exists in Excel",
        mask,
        trans,
        )

    trans = trans.loc[mask].copy()

    # ------------------------------------------------------------------
    # 7. Coordinates
    #
    # LineString:
    #     Use first and last coordinate.
    #
    # Everything else:
    #     Export for manual handling.
    # ------------------------------------------------------------------

    def valid_linestring_coordinates(c):
        return (
            isinstance(c, list)
            and len(c) >= 2
            and isinstance(c[0], (list, tuple))
            and isinstance(c[-1], (list, tuple))
            and len(c[0]) >= 2
            and len(c[-1]) >= 2
        )

    line_mask = trans["geometry_type"].eq("LineString")

    valid_coordinates_mask = (
        line_mask
        & trans["coordinates"].apply(valid_linestring_coordinates)
    )

    print_test_result(
        "usable LineString geometry",
        valid_coordinates_mask,
        trans,
    )

    # Everything that failed the test goes to manual handling
    manual_projects = trans.loc[~valid_coordinates_mask].copy()

    manual_file = OUTPUT_DIR / "manual_projects.csv"

    manual_projects.to_csv(
        manual_file,
        index=False,
    )

    print(
        f"\nManual projects saved to:\n"
        f"  {manual_file}"
    )

    # Keep only automatically processable projects
    trans = trans.loc[valid_coordinates_mask].copy()

    # Extract first and last coordinates
    trans["x0"] = trans["coordinates"].apply(
        lambda c: c[0][0]
    )

    trans["y0"] = trans["coordinates"].apply(
        lambda c: c[0][1]
    )

    trans["x1"] = trans["coordinates"].apply(
        lambda c: c[-1][0]
    )

    trans["y1"] = trans["coordinates"].apply(
        lambda c: c[-1][1]
    )

    # ------------------------------------------------------------------
    # 8. Split AC and DC
    # ------------------------------------------------------------------

    ac_mask = trans["Technology"].eq("AC")
    dc_mask = trans["Technology"].eq("DC")

    print(f"\nAC investments: {ac_mask.sum()}")

    print(f"DC investments: {dc_mask.sum()}")

    print(f"Other technologies: {(~(ac_mask | dc_mask)).sum()}")

    return trans.loc[ac_mask].copy(), trans.loc[dc_mask].copy()

def create_lines_dataframe(trans):
    """Create new_lines.csv dataframe."""

    lines = pd.DataFrame()

    lines["project_status"] = trans.apply(
        create_project_status,
        axis=1,
    )

    lines["length"] = trans["Total route length (km)"]

    lines["build_year"] = (
        trans["Commissioning Year"]
        .astype(str)
        .str.extract(r"(\d{4})")[0]
        .astype("Int64")
    )

    lines["underground"] = False

    lines["v_nom"] = trans["Nominal voltage (kV) "]

    lines["tags"] = trans["Investment number"].apply(create_project_tags)

    lines["x0"] = trans["x0"]
    lines["y0"] = trans["y0"]
    lines["x1"] = trans["x1"]
    lines["y1"] = trans["y1"]

    lines["bus0"] = trans["Substation From"]
    lines["bus1"] = trans["Substation To"]

    # TYNDP2024 project ID
    lines.index = (
        "TYNDP2024_"
        + trans["Investment number"].astype(str)
    )

    lines.index.name = None

    return lines

def create_links_dataframe(trans):
    """Create new_links.csv dataframe."""

    links = pd.DataFrame()

    links["project_status"] = trans.apply(
        create_project_status,
        axis=1,
    )

    links["length"] = trans["Total route length (km)"]

    links["build_year"] = (
        trans["Commissioning Year"]
        .astype(str)
        .str.extract(r"(\d{4})")[0]
        .astype("Int64")
    )

    links["underground"] = False

    links["p_nom"] = trans["Capacity of the investment (MW) "]

    links["tags"] = trans["Investment number"].apply(create_project_tags)

    links["x0"] = trans["x0"]
    links["y0"] = trans["y0"]
    links["x1"] = trans["x1"]
    links["y1"] = trans["y1"]

    links["bus0"] = trans["Substation From"]
    links["bus1"] = trans["Substation To"]

    # TYNDP2024 project ID
    links.index = (
        "TYNDP2024_"
        + trans["Investment number"].astype(str)
    )

    links.index.name = None

    return links

def save_projects(lines, links, output_dir):
    """Save new_lines.csv and new_links.csv."""

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    lines_file = output_dir / "new_lines.csv"
    links_file = output_dir / "new_links.csv"

    lines.to_csv(
        lines_file,
        index=True,
    )

    links.to_csv(
        links_file,
        index=True,
    )

    print("\n" + "=" * 80)
    print("FILES CREATED")
    print("=" * 80)

    print(f"\n{lines_file}")
    print(f"  {len(lines)} lines")

    print(f"\n{links_file}")
    print(f"  {len(links)} links")

def create_project_tags(investment_number):
    """
    Create PyPSA-Eur-style tags for a TYNDP 2024 investment.
    """

    investment_number = int(investment_number)

    return (
        f'"""url""=>""https://tyndp2024.entsoe.eu/projects-map/'
        f'transmission/{investment_number}"", '
        f'""tyndp2024_invest_id""=>""{investment_number}"""'
    )

def create_project_status(row):
    """
    Convert the TYNDP 2024 numeric Status ID into the corresponding
    PyPSA-Eur project status string.

    The mapping is defined in the TYNDP 2024 Status ID column name.
    """

    status_column = next(
        col for col in row.index
        if col.startswith("Status ID")
    )

    status_id = row[status_column]

    status_map = {
        1: "under_consideration",
        2: "in_planning",
        3: "in_permitting",
        4: "under_construction",
        5: "commissioned",
        6: "completed",
    }

    return status_map.get(status_id)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    print("\n" + "=" * 70)
    print("  TYNDP 2024 → PyPSA-Eur Transmission Projects")
    print("=" * 70)

    while True:
        try:
            max_build_year = int(
                input(
                    "\nEnter maximum commissioning year "
                    "(e.g. 2030): "
                )
            )

            if max_build_year < 2000:
                print("Please enter a valid year (>= 2000).")
                continue

            break

        except ValueError:
            print("Please enter the year as a number, e.g. 2030.")

    print(
        f"\n→ Selecting transmission projects commissioned "
        f"by {max_build_year}."
    )

    if not OUTPUT_DIR.exists():
        OUTPUT_DIR.mkdir(
            parents=True,
            exist_ok=True,
        )

    print(f"Created output directory:\n  {OUTPUT_DIR}")

    investments = read_investments(EXCEL_FILE)

    map_data = download_map()

    map_df = create_map_dataframe(map_data)

    trans_lines, trans_links = create_projects(investments, map_df, max_build_year)

    # --------------------------------------------------------------
    # Convert to PyPSA-Eur format
    # --------------------------------------------------------------

    lines = create_lines_dataframe(trans_lines)

    links = create_links_dataframe(trans_links)

    # --------------------------------------------------------------
    # Save
    # --------------------------------------------------------------

    save_projects(
        lines,
        links,
        OUTPUT_DIR,
    )


if __name__ == "__main__":
    main()
