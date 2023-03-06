# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import contextlib
import logging
import os
import sys
from pathlib import Path

import pandas as pd
import pytz
import yaml
from pypsa.components import component_attrs, components
from pypsa.descriptors import Dict
from snakemake.utils import update_config

logger = logging.getLogger(__name__)


# Define a context manager to temporarily mute print statements
@contextlib.contextmanager
def mute_print():
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull):
            yield


def override_component_attrs(directory):
    """Tell PyPSA that links can have multiple outputs by
    overriding the component_attrs. This can be done for
    as many buses as you need with format busi for i = 2,3,4,5,....
    See https://pypsa.org/doc/components.html#link-with-multiple-outputs-or-inputs

    Parameters
    ----------
    directory : string
        Folder where component attributes to override are stored
        analogous to ``pypsa/component_attrs``, e.g. `links.csv`.

    Returns
    -------
    Dictionary of overridden component attributes.
    """

    attrs = Dict({k: v.copy() for k, v in component_attrs.items()})

    for component, list_name in components.list_name.items():
        fn = f"{directory}/{list_name}.csv"
        if os.path.isfile(fn):
            overrides = pd.read_csv(fn, index_col=0, na_values="n/a")
            attrs[component] = overrides.combine_first(attrs[component])

    return attrs


# from pypsa-eur/_helpers.py
def progress_retrieve(url, file):
    import urllib

    from progressbar import ProgressBar

    pbar = ProgressBar(0, 100)

    def dlProgress(count, blockSize, totalSize):
        pbar.update(int(count * blockSize * 100 / totalSize))

    urllib.request.urlretrieve(url, file, reporthook=dlProgress)


def generate_periodic_profiles(dt_index, nodes, weekly_profile, localize=None):
    """
    Give a 24*7 long list of weekly hourly profiles, generate this for each
    country for the period dt_index, taking account of time zones and summer
    time.
    """

    weekly_profile = pd.Series(weekly_profile, range(24 * 7))

    week_df = pd.DataFrame(index=dt_index, columns=nodes)

    for node in nodes:
        timezone = pytz.timezone(pytz.country_timezones[node[:2]][0])
        tz_dt_index = dt_index.tz_convert(timezone)
        week_df[node] = [24 * dt.weekday() + dt.hour for dt in tz_dt_index]
        week_df[node] = week_df[node].map(weekly_profile)

    week_df = week_df.tz_localize(localize)

    return week_df


def parse(l):
    if len(l) == 1:
        return yaml.safe_load(l[0])
    else:
        return {l.pop(0): parse(l)}


def update_config_with_sector_opts(config, sector_opts):
    for o in sector_opts.split("-"):
        if o.startswith("CF+"):
            l = o.split("+")[1:]
            update_config(config, parse(l))
