from math import radians, sin, cos, sqrt, atan2
from pathlib import Path
import pandas as pd



def coordinate_distance_km(lon1, lat1, lon2, lat2):
    """Calculate approximate distance between two lon/lat points in km."""
    earth_radius_km = 6371.0

    lon1, lat1, lon2, lat2 = map(
        radians, [lon1, lat1, lon2, lat2]
        )

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = (
        sin(dlat / 2) ** 2
        + cos(lat1)
        * cos(lat2)
        * sin(dlon / 2) ** 2
        )

    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    return earth_radius_km * c

def create_filtered_transmission_projects(
    base_dir="data/transmission_projects",
    output_folder="personally_filtered_projects",
    max_build_year=2030,
    coordinate_tolerance_km=1.0
    ):

    base_dir = Path(base_dir)
    output_dir = base_dir / output_folder

    file_types = [
        "new_lines.csv",
        "new_links.csv",
        "upgraded_lines.csv",
        "upgraded_links.csv"
        ]

    # Don't include output/template folders as input
    source_dirs = [
        folder
        for folder in base_dir.iterdir()
        if folder.is_dir() and folder.name not in {output_folder, "template", "tyndp2020"}
        ]

    output_dir.mkdir(exist_ok=True)

    print("Source folders:")
    for folder in source_dirs:
        print(f"  {folder}")
    print()

    for filename in file_types:

        dfs = []

        # ---------------------------------------------------------
        # Read files from all source folders
        # ---------------------------------------------------------
        for folder in source_dirs:

            filepath = folder / filename

            if not filepath.exists():
                continue

            # Read original first column as index
            df = pd.read_csv(filepath, index_col=0)

            if "build_year" not in df.columns:
                raise ValueError(f"'build_year' not found in {filepath}")

            # Make sure required coordinate columns exist
            required_coordinates = ["x0", "y0", "x1", "y1"]

            missing = [col for col in required_coordinates if col not in df.columns]

            if missing:
                raise ValueError(f"Missing coordinate columns {missing} in {filepath}")

            # Convert relevant columns to numeric
            df["build_year"] = pd.to_numeric(df["build_year"], errors="coerce")

            for col in required_coordinates:
                df[col] = pd.to_numeric(df[col], errors="coerce")

            # Filter by build year
            df = df[df["build_year"] <= max_build_year].copy()

            if not df.empty:
                dfs.append(df)

            print(
                f"{filename} | {folder.name}: "
                f"{len(df)} projects with "
                f"build_year <= {max_build_year}"
                )

        if not dfs:
            print(f"\n{filename}: no projects found\n")
            continue

        # ---------------------------------------------------------
        # Merge all source folders
        # ---------------------------------------------------------
        merged = pd.concat(dfs, axis=0)

        n_before = len(merged)

        # ---------------------------------------------------------
        # Check duplicate project IDs
        # ---------------------------------------------------------
        duplicate_indices = merged.index[merged.index.duplicated(keep=False)].unique()

        if len(duplicate_indices) > 0:

            print(f"\n{filename}: {len(duplicate_indices)} duplicate project IDs found.")

            for project_id in duplicate_indices:

                duplicate_rows = merged.loc[[project_id]]

                # If the complete rows are identical, safe duplicate
                if duplicate_rows.drop_duplicates().shape[0] == 1:
                    continue

                print("\nCONFLICTING DUPLICATE ID:")
                print(f"Project ID: {project_id}")
                print(duplicate_rows)

                raise ValueError(
                    f"Project ID '{project_id}' occurs multiple times "
                    f"with different data in {filename}."
                    )

            # Remove identical duplicate IDs
            merged = merged[~merged.index.duplicated(keep="first")]

            print("Identical duplicate project IDs removed.")

        # ---------------------------------------------------------
        # Check duplicate projects based on coordinates
        # ---------------------------------------------------------
        coordinate_duplicates = []

        indices = merged.index.tolist()

        for i, idx in enumerate(indices):

            row_i = merged.loc[idx]

            for j in range(i + 1, len(indices)):

                row_j = merged.loc[indices[j]]

                # Skip if coordinates are missing
                if (
                    pd.isna(row_i["x0"])
                    or pd.isna(row_i["y0"])
                    or pd.isna(row_i["x1"])
                    or pd.isna(row_i["y1"])
                    or pd.isna(row_j["x0"])
                    or pd.isna(row_j["y0"])
                    or pd.isna(row_j["x1"])
                    or pd.isna(row_j["y1"])
                    ):
                    continue

                # Same direction
                same_direction = (
                    coordinate_distance_km(
                        row_i["x0"], row_i["y0"],
                        row_j["x0"], row_j["y0"]
                        ) <= coordinate_tolerance_km
                    and
                    coordinate_distance_km(
                        row_i["x1"], row_i["y1"],
                        row_j["x1"], row_j["y1"]
                        ) <= coordinate_tolerance_km
                    )

                # Reversed direction
                reversed_direction = (
                    coordinate_distance_km(
                        row_i["x0"], row_i["y0"],
                        row_j["x1"], row_j["y1"]
                    ) <= coordinate_tolerance_km
                    and
                    coordinate_distance_km(
                        row_i["x1"], row_i["y1"],
                        row_j["x0"], row_j["y0"]
                        ) <= coordinate_tolerance_km
                    )

                if same_direction or reversed_direction:
                    coordinate_duplicates.append((indices[i], indices[j]))

        # ---------------------------------------------------------
        # Report coordinate duplicates
        # ---------------------------------------------------------
        if coordinate_duplicates:

            print(
                f"\n{filename}: "
                f"{len(coordinate_duplicates)} "
                f"coordinate duplicate pairs found:"
            )

            for id1, id2 in coordinate_duplicates:
                print(f"  {id1} <-> {id2}")

            print(
                "\nThese are NOT automatically removed. "
                "Please inspect them first."
            )

        else:
            print(
                f"{filename}: no coordinate duplicates found."
            )

        # ---------------------------------------------------------
        # Output
        # ---------------------------------------------------------
        n_after = len(merged)

        output_path = output_dir / filename

        merged.to_csv(
            output_path,
            index=True,
            index_label=""
        )

        print(
            f"{filename}: "
            f"{n_before} projects before duplicate-ID removal -> "
            f"{n_after} projects written to {output_path}\n"
        )


if __name__ == "__main__":
    create_filtered_transmission_projects(
        max_build_year=2030
    )
