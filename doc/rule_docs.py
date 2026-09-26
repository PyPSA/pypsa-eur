# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Parse Snakemake rule definitions from `rules/*.smk` for the documentation.

The parser is deliberately textual: it reads the rule blocks without loading
the workflow, so the documentation build needs neither Snakemake nor a
configuration. Expressions are reduced to what a reader needs: file paths
without the `{run}` prefix, configuration keys instead of `config_provider`
calls, and "depends on configuration" where an input is chosen at runtime.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

RULES_DIR = Path(__file__).resolve().parents[1] / "rules"

_RULE_RE = re.compile(r"^(\s*)rule (\w+):\s*$")
_SECTION_RE = re.compile(r"^(\w+):\s*(.*)$")
_RUNTIME_MARKERS = ("lambda", "unpack(", "(w)", "dataset_version(", "[", "if ")
_PATH_PATTERNS = [
    (r'resources\(\s*f?"([^"]+)"', "resources/{}"),
    (r'RESULTS\s*\+\s*f?"([^"]+)"', "results/{}"),
    (r'(?<![\w/(])f?"((?:data|config|cutouts)/[^"]+)"', "{}"),
    (r"rules\.(\w+)\.output", "output of `{}`"),
]


@dataclass
class Rule:
    name: str
    file: str
    summary: str = ""
    script: str | None = None
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    settings: list[str] = field(default_factory=list)
    wildcards: list[str] = field(default_factory=list)

    @property
    def module(self) -> str | None:
        """Import path of the script for mkdocstrings, e.g. `build_cop_profiles.run`."""
        return self.script and self.script.removesuffix(".py").replace("/", ".")


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip())


def _split_entries(text: str) -> list[str]:
    """Split a Snakemake section body at top-level commas, dropping keyword names."""
    entries, depth, quote, start = [], 0, None, 0
    for i, ch in enumerate(text):
        if quote:
            quote = None if ch == quote and text[i - 1] != "\\" else quote
        elif ch in "\"'":
            quote = ch
        elif ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            entries.append(text[start:i])
            start = i + 1
    entries.append(text[start:])
    return [re.sub(r"^\w+\s*=", "", e.strip()).strip() for e in entries if e.strip()]


def _paths(expr: str) -> list[str]:
    """Reduce an input or output expression to displayable file paths."""
    found = [
        re.sub(r"\{[^{}]*[\[.('][^{}]*\}", "{...}", fmt.format(m))
        for pattern, fmt in _PATH_PATTERNS
        for m in re.findall(pattern, expr)
    ]
    if found:
        return found
    if any(marker in expr for marker in _RUNTIME_MARKERS):
        return ["depends on configuration"]
    return re.findall(r'f?"([^"]*)"', expr) or [expr]


def _sections(body: list[str]) -> dict[str, str]:
    """Map each Snakemake keyword of a rule body to its joined content."""
    base = _indent(body[0])
    sections: dict[str, str] = {}
    current = None
    for line in body:
        m = _SECTION_RE.match(line.strip())
        if m and _indent(line) == base:
            current = m.group(1)
            sections[current] = m.group(2)
        elif current and _indent(line) > base:
            sections[current] += " " + line.strip()
    return sections


def _parse_body(rule: Rule, body: list[str]) -> None:
    if body[0].strip().startswith('"""'):
        rule.summary = body[0].strip().strip('"').strip()
    sections = _sections(body)
    for entry in _split_entries(sections.get("input", "")):
        rule.inputs += _paths(entry)
    for entry in _split_entries(sections.get("output", "")):
        rule.outputs += _paths(entry)
    for call in re.findall(r"config_provider\(([^)]*)\)", sections.get("params", "")):
        rule.settings.append(".".join(re.findall(r'"([^"]+)"', call)))
    if m := re.search(r'scripts\(\s*"([^"]+)"', sections.get("script", "")):
        rule.script = m.group(1)
    rule.inputs = list(dict.fromkeys(rule.inputs))
    rule.outputs = list(dict.fromkeys(rule.outputs))
    rule.settings = list(dict.fromkeys(rule.settings))
    found = {w for p in rule.outputs for w in re.findall(r"\{(\w+)\}", p)}
    rule.wildcards = sorted(found - {"run"})


def parse_file(path: Path) -> list[Rule]:
    lines = path.read_text().split("\n")
    rules: list[Rule] = []
    i = 0
    while i < len(lines):
        if not (m := _RULE_RE.match(lines[i])):
            i += 1
            continue
        indent = len(m.group(1))
        j = i + 1
        while j < len(lines) and (not lines[j].strip() or _indent(lines[j]) > indent):
            j += 1
        rule = Rule(name=m.group(2), file=path.name)
        body = [
            line
            for line in lines[i + 1 : j]
            if line.strip() and not line.lstrip().startswith("#")
        ]
        _parse_body(rule, body)
        rules.append(rule)
        i = j
    return rules


def parse_rules(files: list[str] | None = None) -> dict[str, Rule]:
    """Return all rules keyed by name; a name defined twice keeps the first definition."""
    paths = [RULES_DIR / f for f in files] if files else sorted(RULES_DIR.glob("*.smk"))
    rules: dict[str, Rule] = {}
    for path in paths:
        for rule in parse_file(path):
            rules.setdefault(rule.name, rule)
    return rules
