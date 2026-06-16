# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import json
import re
from functools import lru_cache
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.default.yaml"
PLOTTING_CONFIG_PATH = ROOT / "config" / "plotting.default.yaml"
SCHEMA_PATH = ROOT / "config" / "schema.default.json"


def _dump(node):
    text = yaml.safe_dump(node, sort_keys=False, default_flow_style=False, width=4096)
    return text.rstrip("\n")


def _schema_type(prop):
    if "enum" in prop:
        return "enum (" + ", ".join(f"`{e}`" for e in prop["enum"]) + ")"
    if "anyOf" in prop:
        return " \\| ".join(_schema_type(b) for b in prop["anyOf"])
    t = prop.get("type")
    if t == "array":
        return "list of " + _schema_type(prop.get("items", {}))
    if "properties" in prop:
        return "any"
    if t == "object":
        ap = prop.get("additionalProperties")
        return f"dict (str -> {_schema_type(ap)})" if isinstance(ap, dict) else "object"
    if isinstance(t, list):
        return " \\| ".join(t)
    return t or "any"


def _schema_default(prop):
    d = prop.get("default")
    if "default" not in prop or d is None:
        return ""
    if isinstance(d, bool):
        return "`true`" if d else "`false`"
    return '`""`' if d == "" else f"`{d}`"


def _schema_desc(prop):
    # collapse RST inline literals/roles left in plain descriptions to markdown
    d = prop.get("markdownDescription", prop.get("description", ""))
    return re.sub(r":math:`([^`]+)`|``([^`]+)``", lambda m: f"`{m[1] or m[2]}`", d)


def _schema_children(prop):
    branches = prop.get("anyOf", [prop])
    return next((b["properties"] for b in branches if "properties" in b), None)


def _schema_rows(props, depth=0):
    prefix = "↳" * depth + (" " if depth else "")
    rows = []
    for key, prop in props.items():
        rows.append(
            f"| {prefix}`{key}` | {_schema_type(prop)} | "
            f"{_schema_default(prop)} | {_schema_desc(prop)} |"
        )
        children = _schema_children(prop)
        if children:
            rows += _schema_rows(children, depth + 1)
    return rows


@lru_cache(maxsize=1)
def _schema():
    return json.loads(SCHEMA_PATH.read_text())


def _resolve(data, path):
    node = data
    for part in path:
        node = node[part]
    return node


def define_env(env):
    _cache = {}

    def _load(source):
        if source not in _cache:
            if source == "plotting":
                path = PLOTTING_CONFIG_PATH
            elif source == "config":
                path = CONFIG_PATH
            else:
                path = ROOT / "config" / source
            _cache[source] = yaml.safe_load(path.read_text())
        return _cache[source]

    @env.macro
    def schema_table(path):
        """Render the config schema at ``path`` as a Property/Type/Default/Description table."""
        node = _schema()
        for part in path.split("."):
            node = node["properties"][part]
        props = node.get("properties", {})
        if not props:
            raise ValueError(f"schema path '{path}' has no properties to tabulate")
        header = (
            "| Property | Type | Default | Description |\n"
            "|----------|------|---------|-------------|"
        )
        return header + "\n" + "\n".join(_schema_rows(props))

    @env.macro
    def yaml_section(*paths, source="config", with_key=True):
        data = _load(source)

        if not paths:
            return _dump(data)

        parts_list = [p.split(".") for p in paths]

        if len(paths) == 1:
            parts = parts_list[0]
            node = _resolve(data, parts)
            if with_key:
                return _dump({parts[-1]: node})
            return _dump(node)

        parents = {tuple(p[:-1]) for p in parts_list}
        if len(parents) != 1:
            raise ValueError("paths must share a common parent")
        parent_path = list(parents.pop())
        children = [p[-1] for p in parts_list]
        parent_node = _resolve(data, parent_path)
        subset = {ck: parent_node[ck] for ck in children}
        if with_key and parent_path:
            return _dump({parent_path[-1]: subset})
        return _dump(subset)
