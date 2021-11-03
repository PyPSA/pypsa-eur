"""
Preprocess gas network based on data from:

    [1] the SciGRID Gas project
        (https://www.gas.scigrid.de/)

    [2] ENTSOG capacity map
        (https://www.entsog.eu/sites/default/files/2019-10/Capacities%20for%20Transmission%20Capacity%20Map%20RTS008_NS%20-%20DWH_final.xlsx)
"""

import logging
logger = logging.getLogger(__name__)

import re
import json

import pandas as pd
from shapely.geometry import Point
from pypsa.geo import haversine_pts


def string2list(string, with_none=True):
    """Convert string format to a list."""

    if with_none:
        p2 = re.compile('None')
        string = p2.sub('\"None\"', string)
    else:
        p = re.compile('(?<!\\\\)\'')
        string = p.sub('\"', string)

    return json.loads(string)


def diameter2capacity(pipe_diameter_mm):
    """Calculate pipe capacity in MW based on diameter in mm.

    20 inch (500 mm)  50 bar -> 1.5   GW CH4 pipe capacity (LHV)
    24 inch (600 mm)  50 bar -> 5     GW CH4 pipe capacity (LHV)
    36 inch (900 mm)  50 bar -> 11.25 GW CH4 pipe capacity (LHV)
    48 inch (1200 mm) 80 bar -> 21.7  GW CH4 pipe capacity (LHV)

    Based on p.15 of https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf
    """

    # slopes definitions
    m0 = (1500 - 0) / (500 - 0)
    m1 = (5000 - 1500) / (600 - 500)
    m2 = (11250 - 5000) / (900 - 600)
    m3 = (21700 - 11250) / (1200 - 900)

    # intercept
    a0 = 0
    a1 = -16000
    a2 = -7500
    a3 = -20100

    if pipe_diameter_mm < 500:
        return a0 + m0 * pipe_diameter_mm
    elif pipe_diameter_mm < 600:
        return a1 + m1 * pipe_diameter_mm
    elif pipe_diameter_mm < 900:
        return a2 + m2 * pipe_diameter_mm
    else:
        return a3 + m3 * pipe_diameter_mm


def find_terminal_points(df):
    
    latlon = []

    for attr in ["lat", "long"]:
    
        s = df[attr].apply(string2list)

        s = s.apply(lambda x: [x[0], x[-1]])

        latlon.append(pd.DataFrame(s.to_list(),
            columns=[f"{attr}0", f"{attr}1"]
        ))
    
    latlon = pd.concat(latlon, axis=1)
    
    points = latlon.apply(
        lambda x: {
            "point0": Point(x.long0, x.lat0),
            "point1": Point(x.long1, x.lat1)
        },
        axis=1,
        result_type='expand'
    )
    
    return pd.concat([df, points], axis=1)


def process_gas_network_data(fn):

    df = pd.read_csv(fn, sep=',')

    df = find_terminal_points(df)

    to_drop = ["name", "source_id", "country_code", "node_id",
               "long", "lat", "lat_mean", "long_mean", "num_compressor"]
    df.drop(to_drop, axis=1, inplace=True)

    to_rename = {
        "is_bothDirection": "bidirectional",
        "start_year": "build_year",
        "length_km": "length",
        "Capacity_GWh_h": "p_nom_data",
        "id": "tags",
    }
    df.rename(columns=to_rename, inplace=True)
    
    df.bidirectional = df.bidirectional.astype(bool)

    # convert from GWh/h to MW
    df.p_nom_data *= 1e3

    # for pipes with missing diameter, assume 500 mm
    df.loc[df.diameter_mm.isna(), "diameter_mm"] = 500.

    # for nord stream and small pipelines take original capacity data
    # otherwise inferred values from pipe diameter
    df["p_nom"] = df.diameter_mm.map(diameter2capacity)
    df.p_nom.update(
        df.p_nom_data.where((df.diameter_mm < 500) | (df.max_pressure_bar == 220))
    )

    df["length_haversine"] = df.apply(
        lambda p: 1.5 * haversine_pts([p.point0.x, p.point1.y], [p.point1.x, p.point1.y]),
        axis=1
    )

    df.length.update(df.length_haversine.where(df.length.isna()))
    
    return df


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_gas_network')

    logging.basicConfig(level=snakemake.config['logging_level'])

    gas_network = process_gas_network_data(snakemake.input.gas_network)

    gas_network.to_csv(snakemake.output.cleaned_gas_network)