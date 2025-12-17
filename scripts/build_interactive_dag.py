# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Generate an interactive HTML visualization of the Snakemake rulegraph
with a right-hand sidebar showing detailed filegraph metadata extracted
from Graphviz HTML labels.
"""

import json
import re
import textwrap
import urllib.request

import pydot
import yaml

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)


def load_linguist_extension_colors() -> dict[str, str]:
    """
    Fetch GitHub Linguist's languages.yml and return extension→color mapping.
    Returns an empty dict on any network error, YAML error, or unexpected type.
    """
    url = (
        "https://raw.githubusercontent.com/github-linguist/"
        "linguist/main/lib/linguist/languages.yml"
    )

    # Network fetch
    try:
        with urllib.request.urlopen(url, timeout=5) as r:
            raw = r.read()
    except Exception:
        return {}

    # Decode
    try:
        data = raw.decode("utf-8")
    except Exception:
        return {}

    # YAML parse
    try:
        languages = yaml.safe_load(data)
    except Exception:
        return {}

    if not isinstance(languages, dict):
        return {}

    ext_to_color: dict[str, str] = {}

    for attrs in languages.values():
        if not isinstance(attrs, dict):
            continue

        color = attrs.get("color")
        exts = attrs.get("extensions")

        if not isinstance(color, str):
            continue
        if not isinstance(exts, list):
            continue

        for ext in exts:
            if isinstance(ext, str) and ext.startswith("."):
                ext_to_color[ext.lower()] = color

    return ext_to_color


HOVER_DELAY_MS = 20
DEFAULT_INFO_TEXT = "<em>Select a rule to see inputs/outputs.</em>"
COLOR_ACTIVE_NODE = "#d62728"
COLOR_UPSTREAM_NODE = "#ff7f0e"
COLOR_DOWNSTREAM_NODE = "#1f77b4"
COLOR_EDGE_HOVER = "#2ca02c"
COLOR_BACKGROUND = "#ffffff"
COLOR_SIDEBAR_BG = "#f9f9f9"
COLOR_SIDEBAR_BORDER = "#ddd"
COLOR_WILDCARD = "#1f77ff"
FILE_EXTENSION_COLORS = load_linguist_extension_colors()
COLOR_EXT_DEFAULT = "#aaaaaa"


def clean_svg(svg_text: str) -> str:
    cleaned = []
    skip = False
    for line in svg_text.splitlines():
        s = line.strip()
        if s.startswith("<?xml"):
            continue
        if s.startswith("<!DOCTYPE"):
            skip = True
            continue
        if skip:
            if s.endswith(">"):
                skip = False
            continue
        cleaned.append(line)

    return "\n".join(cleaned)


def _highlight_extensions(html: str) -> str:
    if not html:
        return html

    def color_last_ext(filename: str) -> str:
        before, _, tail = filename.rpartition("/")
        last = tail

        # No extension → nothing to color
        if "." not in last:
            return filename

        # Identify last dot
        dot_index = last.rfind(".")
        ext = last[dot_index:].lower()  # e.g. ".xlsx"

        # Determine color: Linguist if known, fallback otherwise
        color = FILE_EXTENSION_COLORS.get(ext, COLOR_EXT_DEFAULT)

        colored_last = (
            last[:dot_index]
            + f'<span class="ext" style="color:{color}; font-weight:bold;">{last[dot_index:]}</span>'
        )

        return before + ("/" if before else "") + colored_last

    # Split HTML into tags "<...>" and text nodes
    parts = re.split(r"(<[^>]+>)", html)

    out = []
    for part in parts:
        if part.startswith("<") and part.endswith(">"):
            out.append(part)
        else:
            # color filename-like tokens in text nodes
            tokens = part.split()
            new_tokens = []
            for tok in tokens:
                if "." in tok:
                    new_tokens.append(color_last_ext(tok))
                else:
                    new_tokens.append(tok)
            out.append(" ".join(new_tokens))

    return "".join(out)


def _highlight_wildcards(html: str) -> str:
    if not html:
        return html

    out = []
    i = 0
    n = len(html)

    while i < n:
        if html[i] == "{":
            start = i
            depth = 1
            i += 1
            while i < n and depth > 0:
                if html[i] == "{":
                    depth += 1
                elif html[i] == "}":
                    depth -= 1
                i += 1
            segment = html[start:i]
            out.append(f'<span class="wc">{segment}</span>')
        else:
            out.append(html[i])
            i += 1

    return "".join(out)


def style_html_table(html: str) -> str:
    if not html:
        return html

    html = re.sub(r"<hr\s*/?>", "", html, flags=re.IGNORECASE)
    html = _highlight_extensions(html)
    html = _highlight_wildcards(html)
    return html


def parse_filegraph_html_tables(graph: pydot.Dot) -> dict[str, str]:
    def extract_html(raw: str | None) -> str | None:
        if not raw:
            return None
        s = raw.strip()
        if s.startswith("<") and s.endswith(">"):
            s = s[1:-1].strip()
        return s or None

    out = {}
    for node in graph.get_nodes():
        name = node.get_name().strip('"')
        if name in {"node", "edge", "graph"}:
            continue

        raw = node.get_attributes().get("label")
        html = extract_html(raw)
        if html:
            out[name] = style_html_table(html)

    return out


def build_html(svg_embed: str, node_tables: dict[str, str], svg_panzoom_js: str) -> str:
    safe_panzoom_js = svg_panzoom_js.replace("</script>", "<\\/script>")
    tables_json = json.dumps(node_tables)

    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Interactive Rulegraph</title>

        <style>
            html, body {{
                margin: 0;
                padding: 0;
                background: {COLOR_BACKGROUND};
                overflow: hidden;
                font-family: sans-serif;
            }}

            #svg-container {{
                position: absolute;
                left: 0;
                top: 0;
                width: calc(100vw - 340px);
                height: 100vh;
                overflow: hidden;
                background: {COLOR_BACKGROUND};
            }}

            /* Resizable right panel */
            #info-panel {{
                position: absolute;
                right: 0;
                top: 0;
                width: 340px;
                height: 100vh;
                background: {COLOR_SIDEBAR_BG};
                border-left: 1px solid {COLOR_SIDEBAR_BORDER};
                overflow-y: auto;
                overflow-x: hidden;
                padding: 16px;
                box-sizing: border-box;
                resize: horizontal;
                min-width: 250px;
                max-width: 700px;
            }}

            #info-panel table {{
                width: 100%;
                border-collapse: collapse;
                font-size: 13px;
                table-layout: fixed;
                background: white;
            }}

            #info-panel td {{
                padding: 6px 8px;
                border-bottom: 1px solid #e5e5e5;
                word-break: break-word;
            }}

            #info-panel tr:nth-child(even) td {{
                background: #fafafa;
            }}

            .wc {{
                color: {COLOR_WILDCARD};
                font-weight: bold;
            }}

            .selected-node * {{
                fill-opacity: 0.35 !important;
            }}

            svg {{
                background: {COLOR_BACKGROUND};
            }}

            .active-node polygon,
            .active-node ellipse {{
                stroke: {COLOR_ACTIVE_NODE};
                stroke-width: 3px;
            }}

            .upstream-node polygon,
            .upstream-node ellipse {{
                stroke: {COLOR_UPSTREAM_NODE};
                stroke-width: 2.5px;
            }}

            .downstream-node polygon,
            .downstream-node ellipse {{
                stroke: {COLOR_DOWNSTREAM_NODE};
                stroke-width: 2.5px;
            }}

            .edge-hover path {{
                stroke: {COLOR_EDGE_HOVER};
                stroke-width: 3px;
            }}

            .upstream-edge path {{
                stroke: {COLOR_UPSTREAM_NODE};
                stroke-width: 2.5px;
            }}

            .downstream-edge path {{
                stroke: {COLOR_DOWNSTREAM_NODE};
                stroke-width: 2.5px;
            }}

            .dimmed * {{
                opacity: 0.1;
            }}
            .isolated-only * {{
                opacity: 1;
            }}

            .node-hover * {{
                fill-opacity: 0.25;
            }}
        </style>
    </head>

    <body>

        <div id="svg-container">
            {svg_embed}
        </div>

        <div id="info-panel">
            <div id="info-panel-content">{DEFAULT_INFO_TEXT}</div>
        </div>

        <script>const nodeTables = {tables_json};</script>

        <script>
            {safe_panzoom_js}
        </script>

        <script>
        (function() {{
            const sidebar = document.getElementById("info-panel-content");
            const container = document.getElementById("svg-container");
            const svg = container.querySelector("svg");
            if (!svg) return;

            const DEFAULT_INFO = `{DEFAULT_INFO_TEXT}`;
            const HOVER_DELAY = {HOVER_DELAY_MS};
            let hoverTimer = null;
            let isolated = null;
            let selectedNode = null;

            const originalFill = new WeakMap();

            const panZoom = svgPanZoom(svg, {{
                controlIconsEnabled: false,
                minZoom: 0.1,
                maxZoom: 20,
                zoomScaleSensitivity: 0.3,
                fit: true,
                center: true
            }});

            const nodes = Array.from(svg.querySelectorAll("g.node"));
            const edges = Array.from(svg.querySelectorAll("g.edge"));

            const nodeEl = new Map();
            const upstream = new Map();
            const downstream = new Map();

            // Build node registry
            nodes.forEach(g => {{
                const t = g.querySelector("title");
                if (!t) return;
                const id = t.textContent.trim();
                nodeEl.set(id, g);
                upstream.set(id, new Set());
                downstream.set(id, new Set());

                g.setAttribute("pointer-events", "bounding-box");
                const shape = g.querySelector("polygon, ellipse, path");
                if (shape) {{
                    shape.setAttribute("pointer-events", "all");
                    if (!shape.getAttribute("fill") || shape.getAttribute("fill") === "none")
                        shape.setAttribute("fill", "transparent");
                }}
            }});

            // Build adjacency
            edges.forEach(edge => {{
                const t = edge.querySelector("title");
                if (!t) return;
                const [src, dst] = t.textContent.trim().split("->").map(s => s.trim());
                if (nodeEl.has(src) && nodeEl.has(dst)) {{
                    downstream.get(src).add(dst);
                    upstream.get(dst).add(src);
                }}
            }});

            function bfs(start, adj) {{
                const visited = new Set();
                const stack = [start];
                while (stack.length) {{
                    const u = stack.pop();
                    for (const v of (adj.get(u) || [])) {{
                        if (!visited.has(v)) {{
                            visited.add(v);
                            stack.push(v);
                        }}
                    }}
                }}
                return visited;
            }}

            function reset_selection() {{
                if (selectedNode) {{
                    selectedNode.classList.remove("selected-node");
                }}
                selectedNode = null;
            }}

            function clearAll() {{
                svg.querySelectorAll(
                    ".active-node,.upstream-node,.downstream-node," +
                    ".upstream-edge,.downstream-edge,.edge-hover," +
                    ".dimmed,.isolated-only,.node-hover"
                ).forEach(el => el.classList.remove(
                    "active-node","upstream-node","downstream-node",
                    "upstream-edge","downstream-edge","edge-hover",
                    "dimmed","isolated-only","node-hover"
                ));

                if (!selectedNode) {{
                    nodes.forEach(g => {{
                        if (g.classList.contains("no-restore-fill")) return;
                        const shape = g.querySelector("polygon, ellipse, path");
                        if (shape && originalFill.has(shape)) {{
                            shape.setAttribute("fill", originalFill.get(shape));
                        }}
                    }});
                }}
            }}

            function isolateNode(id) {{
                clearAll();
                isolated = id;

                const ups = bfs(id, upstream);
                const downs = bfs(id, downstream);
                const keep = new Set([id, ...ups, ...downs]);

                svg.querySelectorAll("g.node, g.edge").forEach(el => el.classList.add("dimmed"));
                keep.forEach(k => nodeEl.get(k)?.classList.add("isolated-only"));

                edges.forEach(edge => {{
                    const t = edge.querySelector("title");
                    if (!t) return;
                    const [src, dst] = t.textContent.trim().split("->").map(s => s.trim());
                    if (keep.has(src) && keep.has(dst))
                        edge.classList.add("isolated-only");
                }});
            }}

            function highlight(id) {{
                clearAll();

                const ups = bfs(id, upstream);
                const downs = bfs(id, downstream);

                nodeEl.get(id)?.classList.add("active-node");
                ups.forEach(u => nodeEl.get(u)?.classList.add("upstream-node"));
                downs.forEach(d => nodeEl.get(d)?.classList.add("downstream-node"));

                edges.forEach(edge => {{
                    const t = edge.querySelector("title");
                    if (!t) return;
                    const [src, dst] = t.textContent.trim().split("->").map(s => s.trim());
                    if (downs.has(dst) && (downs.has(src) || src === id))
                        edge.classList.add("downstream-edge");
                    if (ups.has(src) && (ups.has(dst) || dst === id))
                        edge.classList.add("upstream-edge");
                }});

                svg.querySelectorAll("g.node, g.edge").forEach(el => {{
                    if (!el.classList.contains("active-node") &&
                        !el.classList.contains("upstream-node") &&
                        !el.classList.contains("downstream-node") &&
                        !el.classList.contains("upstream-edge") &&
                        !el.classList.contains("downstream-edge"))
                        el.classList.add("dimmed");
                }});
            }}

            // Edge hover
            edges.forEach(edge => {{
                edge.addEventListener("mouseenter", () => {{
                    clearTimeout(hoverTimer);
                    hoverTimer = setTimeout(() => {{
                        if (!edge.matches(":hover")) return;

                        clearAll();
                        edge.classList.add("edge-hover");

                        const [src, dst] = edge.querySelector("title")
                            .textContent.trim().split("->").map(s => s.trim());

                        nodeEl.get(src)?.classList.add("active-node");
                        nodeEl.get(dst)?.classList.add("active-node");

                        svg.querySelectorAll("g.node, g.edge").forEach(el => {{
                            if (!el.classList.contains("active-node") &&
                                !el.classList.contains("edge-hover"))
                                el.classList.add("dimmed");
                        }});
                    }}, HOVER_DELAY);
                }});

                edge.addEventListener("mouseleave", () => {{
                    clearTimeout(hoverTimer);
                    if (!isolated) clearAll();
                    else isolateNode(isolated);
                }});
            }});

            // Node hover + click
            nodeEl.forEach((g, id) => {{
                const shape = g.querySelector("polygon, ellipse, path");

                g.addEventListener("mouseenter", () => {{
                    clearTimeout(hoverTimer);
                    hoverTimer = setTimeout(() => {{
                        if (!g.matches(":hover")) return;

                        if (!g.classList.contains("selected-node")) {{
                            if (shape) {{
                                if (!originalFill.has(shape))
                                    originalFill.set(shape, shape.getAttribute("fill") || "none");
                                const stroke = shape.getAttribute("stroke") || "#000";
                                shape.setAttribute("fill", stroke);
                                g.classList.add("node-hover");
                                g.classList.add("no-restore-fill");
                            }}
                        }}

                        if (!isolated && !g.classList.contains("selected-node"))
                            highlight(id);
                    }}, HOVER_DELAY);
                }});

                g.addEventListener("mouseleave", () => {{
                    clearTimeout(hoverTimer);

                    if (!g.classList.contains("selected-node")) {{
                        if (shape) {{
                            g.classList.remove("node-hover");
                            g.classList.remove("no-restore-fill");
                            shape.setAttribute("fill", originalFill.get(shape) || "none");
                        }}
                    }}

                    if (!isolated)
                        clearAll();
                    else
                        isolateNode(isolated);
                }});

                g.addEventListener("click", (e) => {{
                    e.stopPropagation();
                    clearTimeout(hoverTimer);

                    reset_selection();

                    selectedNode = g;
                    selectedNode.classList.add("selected-node");

                    if (isolated === id) {{
                        isolated = null;
                        clearAll();
                        sidebar.innerHTML = DEFAULT_INFO;
                    }} else {{
                        isolateNode(id);
                        sidebar.innerHTML = nodeTables[id] || "<em>No metadata available.</em>";
                    }}
                }});
            }});

            // Clicking background clears
            svg.addEventListener("click", () => {{
                isolated = null;
                reset_selection();
                clearAll();
                sidebar.innerHTML = DEFAULT_INFO;
            }});

            // Auto-fit view
            requestAnimationFrame(() => {{
                panZoom.resize();
                panZoom.updateBBox();

                const bbox = svg.getBBox();
                const zoom = Math.min(
                    container.clientWidth / bbox.width,
                    container.clientHeight / bbox.height
                );
                panZoom.zoom(zoom);

                const panX = (container.clientWidth - bbox.width * zoom) / 2 - bbox.x * zoom;
                const panY = (container.clientHeight - bbox.height * zoom) / 2 - bbox.y * zoom;

                panZoom.pan({{x: panX, y: panY}});
            }});

        }})();
        </script>

    </body>
    </html>
    """

    return textwrap.dedent(html)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_interactive_dag")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    svg_raw = open(snakemake.input.svg, encoding="utf-8").read()
    svg_embed = clean_svg(svg_raw)

    filegraph_raw = open(snakemake.input.filegraph, encoding="utf-8").read()
    filegraph = pydot.graph_from_dot_data(filegraph_raw)[0]

    # Extract HTML tables from filegraph nodes
    node_tables = parse_filegraph_html_tables(filegraph)

    # Load SVG panzoom JS
    with open(snakemake.input.js, encoding="utf-8") as f:
        svg_panzoom_js = f.read()

    # Build final HTML
    html_output = build_html(svg_embed, node_tables, svg_panzoom_js)

    # Write output HTML
    with open(snakemake.output.html, "w", encoding="utf-8") as f:
        f.write(html_output)
