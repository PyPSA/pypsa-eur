# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Consistency checks between the workflow and its documentation."""

import re
from pathlib import Path

import pytest

from doc.rule_docs import parse_rules

DOCUMENTED_FILES = [
    "build_electricity.smk",
    "build_sector.smk",
    "compose.smk",
    "solve.smk",
    "postprocess.smk",
]


def _documented_rule_names() -> set[str]:
    names: set[str] = set()
    for page in Path("doc/rules").glob("*.md"):
        text = page.read_text()
        for call in re.findall(r"\{\{\s*rules\((.*?)\)\s*\}\}", text, re.S):
            names.update(re.findall(r'"(\w+)"', call))
    return names


def test_every_rule_is_documented():
    """Each rule of the build, compose, solve and postprocess files is on a rules page."""
    expected = set(parse_rules(DOCUMENTED_FILES)) | {"build_natura_raster"}
    documented = _documented_rule_names()
    assert expected - documented == set()
    assert documented - expected == set()


@pytest.mark.parametrize("rule", parse_rules().values(), ids=lambda r: r.name)
def test_every_rule_has_summary(rule):
    """Each rule carries a one-sentence docstring."""
    assert rule.summary, f"rule {rule.name} in {rule.file} has no docstring"
    assert rule.summary.endswith("."), rule.summary
    assert len(rule.summary) <= 100, rule.summary


def test_rule_scripts_exist():
    for rule in parse_rules().values():
        if rule.script:
            assert Path("scripts", rule.script).exists(), rule.name
