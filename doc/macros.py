# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.default.yaml"
PLOTTING_CONFIG_PATH = ROOT / "config" / "plotting.default.yaml"


def _dump(node):
    text = yaml.safe_dump(node, sort_keys=False, default_flow_style=False, width=4096)
    return text.rstrip("\n")


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
