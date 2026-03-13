# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from datetime import date
from pathlib import Path

import pandas as pd
import pandera.pandas as pa
from natsort import natsort_keygen
from pandera.pandas import Check, Column

VERSIONS_CSV = Path(__file__).parent.parent / "data" / "versions.csv"
VALID_SOURCES = ["primary", "archive", "build"]  # Order defines sort priority

VALID_TAGS = {
    "latest",
    "supported",
    "not-supported",
    "deprecated",
    "might-work",
    "not-tested",
    "broken-link",
}

not_empty = [Check.str_length(min_value=1), Check.str_matches(r"\S")]
valid_tags = Check(lambda s: all(t in VALID_TAGS for t in s.split()), element_wise=True)
url_safe = Check.str_matches(
    r"^[a-z0-9_\-\.]+$",
    error="Version must be URL-safe (only alphanumeric, hyphen, underscore, dot)",
)


def sort_versions(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by dataset (asc), version (desc, natural), source (as predefined)."""
    natsort_key = natsort_keygen(key=str.casefold)
    df = df.copy()
    df["source"] = pd.Categorical(df["source"], categories=VALID_SOURCES, ordered=True)
    df = df.sort_values(
        by=["dataset", "version", "source"],
        key=lambda c: c.map(natsort_key) if c.name != "source" else c,
        ascending=[True, False, True],
    ).reset_index(drop=True)
    df["source"] = df["source"].astype(str)
    return df


is_sorted = Check(lambda df: df.equals(sort_versions(df)), error="Data must be sorted")
archive_has_url = Check(
    lambda df: df.loc[df["source"] == "archive", "url"].str.len().gt(0).all(),
    error="Archive entries must have a URL",
)
one_latest_per_dataset_source = Check(
    lambda df: df[df["tags"].str.contains("latest")]
    .groupby(["dataset", "source"])
    .size()
    .eq(1)
    .all(),
    error="Exactly one 'latest' tag required per dataset/source combination",
)
latest_same_version_across_sources = Check(
    lambda df: df[(df["tags"].str.contains("latest")) & (df["version"] != "unknown")]
    .groupby("dataset")["version"]
    .nunique()
    .le(1)
    .all(),
    error="All 'latest' entries for a dataset must have the same version across sources (excluding 'unknown')",
)
VersionsSchema = pa.DataFrameSchema(
    {
        "dataset": Column(str, not_empty, nullable=False),
        "version": Column(str, not_empty + [url_safe], nullable=False),
        "source": Column(str, Check.isin(VALID_SOURCES)),
        "tags": Column(str, valid_tags, nullable=False),
        "added": Column(
            str,
            Check.str_matches(r"^\d{4}-\d{2}-\d{2}$"),
            nullable=False,
            coerce=True,
            default=date.today().isoformat(),
        ),
        "note": Column(str, nullable=True),
        "url": Column(
            str,
            Check.str_matches(
                r'^(https?://[^:\s<>"|*]*)?$',
                error='URL must start with http(s):// and not contain colons (use %3A), spaces, or Windows-invalid characters (<>"|*).',
            ),
            nullable=True,
        ),
    },
    checks=[
        is_sorted,
        archive_has_url,
        one_latest_per_dataset_source,
        latest_same_version_across_sources,
    ],
    coerce=True,
    strict=True,
    ordered=True,
    unique=["dataset", "version", "source"],
)


def load_versions() -> pd.DataFrame:
    return pd.read_csv(VERSIONS_CSV, dtype=str, keep_default_na=False)


def validate_versions(fix: bool = False) -> pd.DataFrame:
    df = load_versions()
    if fix:
        df = sort_versions(df)
    try:
        df = VersionsSchema.validate(df, lazy=True)
    except pa.errors.SchemaErrors as e:
        msg = f"{e.message}\n\nTry 'pixi run python test/test_data_versions_layer.py' to auto-fix (sorting, defaults, etc.)."
        e.message = msg
        raise
    if fix:
        df.to_csv(VERSIONS_CSV, index=False)
    return df


def test_versions_csv():
    validate_versions()


if __name__ == "__main__":
    validate_versions(fix=True)
    print("versions.csv validated and fixed")
