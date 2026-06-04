# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import os

import requests


def on_page_markdown(markdown, page, config, files, **kwargs):
    import re

    markdown = re.sub(r"(\[\^[^\]]+\]:.*)\n(\[\^)", r"\1\n\n\2", markdown)
    return markdown


def on_post_build(config, **kwargs):
    url = "https://zenodo.org/records/14144752/files/map.html?download=1"
    out = os.path.join(config["site_dir"], "base-network-raw.html")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write(requests.get(url).text)
