# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from pathlib import Path

import pandas as pd
import pandera.pandas as pa
import pytest
import yaml

from scripts._helpers import _load_data_version, load_data_versions
from scripts.lib.validation.config import (
    validate_config,
)
from scripts.lib.validation.config.data import VersionsSchema, sort_versions

VERSIONS_PATHS = [
    Path(__file__).parent.parent / file if not Path(file).is_absolute() else Path(file)
    for file in validate_config({}).data.version_files
]


@pytest.mark.parametrize("file", VERSIONS_PATHS, ids=[str(p) for p in VERSIONS_PATHS])
def test_versions_csv(pytestconfig, file):
    fix = pytestconfig.getoption("fix")
    df = _load_data_version(file, validate=False)
    if fix:
        df = sort_versions(df)
    try:
        df = VersionsSchema.validate(df, lazy=True)
    except pa.errors.SchemaErrors as e:
        if fix:
            msg = f"{e.message}\n\nAttempted to automatically fix the `{file.name}` version file but above issues still persist."
        else:
            msg = f"{e.message}\n\nTry `pixi run -e test pytest test/test_data_versions_layer.py --fix` to attempt to fix the `{file.name}` version file."
        e.message = msg
        raise
    if fix:
        if file.suffix.lower() in {".yaml", ".yml"}:
            yaml_df = df.to_dict(orient="records")
            file.write_text(yaml.safe_dump(yaml_df, sort_keys=False))
        else:
            df.to_csv(file, index=False)


def test_load_data_versions_combined(tmp_path):
    """Test that load_data_versions correctly combines multiple files."""
    file_1 = tmp_path / "versions_1.csv"
    file_2 = tmp_path / "versions_2.yaml"
    file_1.write_text(
        "dataset,version,source,tags,added,note,url\n"
        "dataset1,1.0,primary,latest supported,2024-01-01,,http://example.com/1\n"
        "dataset2,2.0,archive,supported,2024-01-02,,http://example.com/2\n"
    )
    file_2.write_text(
        "- dataset: dataset1\n"
        "  version: 1.1\n"
        "  source: archive\n"
        "  tags: supported\n"
        "  added: 2024-01-03\n"
        "  note: Updated version\n"
        "  url: http://example.com/1.1\n"
        "- dataset: dataset2\n"
        "  version: 2.0\n"
        "  source: archive\n"
        "  tags: supported\n"
        "  added: 2024-01-03\n"
        "  note: Updated version\n"
        "  url: http://example.com/2.1\n"
        "- dataset: dataset3\n"
        "  version: 3.0\n"
        "  source: primary\n"
        "  tags: latest supported\n"
        "  added: 2024-01-04\n"
        "  note: New dataset\n"
        "  url: http://example.com/3\n"
    )
    combined_df = load_data_versions(file_1, file_2)
    expected_df = pd.DataFrame(
        {
            "dataset": ["dataset1", "dataset1", "dataset2", "dataset3"],
            "version": ["1.0", "1.1", "2.0", "3.0"],
            "source": ["primary", "archive", "archive", "primary"],
            "tags": [
                ["latest", "supported"],
                ["supported"],
                ["supported"],
                ["latest", "supported"],
            ],
            "added": ["2024-01-01", "2024-01-03", "2024-01-03", "2024-01-04"],
            "note": ["", "Updated version", "Updated version", "New dataset"],
            "url": [
                "http://example.com/1",
                "http://example.com/1.1",
                "http://example.com/2.1",
                "http://example.com/3",
            ],
            "latest": [True, False, False, True],
            "supported": [True, True, True, True],
        },
        index=[0, 1, 2, 3],
    )
    pd.testing.assert_frame_equal(combined_df, expected_df)
