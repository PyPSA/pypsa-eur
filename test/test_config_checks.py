# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import warnings
from io import StringIO
from unittest.mock import patch

from scripts._helpers import (
    DeprecationConfigWarning,
    InvalidConfigWarning,
    check_deprecated_config,
    check_invalid_config,
)

SAMPLE_DEPRECATIONS = """
- old_entry: "old:key"
  new_entry: "new:key"

- old_entry: "removed:key"
  message: "This key is obsolete and should be deleted"

- old_entry: "example:old_key"
  new_entry: "example:new_key"
  message: "Custom warning message"
"""


def test_config_deprecations():
    test_config = {
        "old": {"key": "legacy_value"},
        "removed": {"key": "dangerous_value"},
        "example": {"old_key": "original_value", "new_key": "existing_value"},
        "unrelated": {"data": "untouched"},
    }

    with warnings.catch_warnings(record=True) as captured_warnings:
        with patch("builtins.open", return_value=StringIO(SAMPLE_DEPRECATIONS)):
            check_deprecated_config(test_config, "dummy_path.yaml")

    # Verify warnings
    assert len(captured_warnings) == 3

    warning_messages = [str(w.message) for w in captured_warnings]

    # Check basic rename warning
    assert any(
        "'old:key' is deprecated. Use 'new:key' instead" in msg
        for msg in warning_messages
    )

    # Check removal warning with custom message
    assert any("obsolete and should be deleted" in msg for msg in warning_messages)

    # Check custom message and conflict warning
    assert any("Custom warning message" in msg for msg in warning_messages)
    assert any(
        "Both keys present - remove deprecated entry" in msg for msg in warning_messages
    )

    # Verify warning types
    assert all(
        isinstance(w.message, DeprecationConfigWarning) for w in captured_warnings
    )

    # Verify config updates
    assert test_config["new"]["key"] == "legacy_value"  # Renamed value
    assert (
        test_config["example"]["new_key"] == "existing_value"
    )  # Existing value preserved
    assert "key" in test_config["removed"]  # Removed key not deleted (just warned)
    assert test_config["unrelated"] == {"data": "untouched"}  # Unrelated data unchanged


def test_config_invalid_entries():
    test_config = {
        "valid_section": {"nested_valid": "ok"},
        "invalid_section": {"bad_key": "value"},
        "clustering": {"invalid_option": "bad"},
    }

    default_config = """
valid_section:
    nested_valid: default
    other_valid: default
clustering:
    temporal:
        resolution: 1
    """

    with warnings.catch_warnings(record=True) as captured_warnings:
        with patch("builtins.open", return_value=StringIO(default_config)):
            check_invalid_config(test_config, "dummy_default.yaml")

    warning_messages = [str(w.message) for w in captured_warnings]

    # Check warning for invalid top-level section
    assert any(
        "Config entry 'invalid_section' is not supported" in msg
        for msg in warning_messages
    )

    # Check warning for invalid nested option
    assert any(
        "Config entry 'clustering:invalid_option' is not supported" in msg
        for msg in warning_messages
    )

    # Verify warning types
    assert all(isinstance(w.message, InvalidConfigWarning) for w in captured_warnings)
