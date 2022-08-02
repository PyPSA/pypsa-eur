

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
plt.style.use('ggplot')

from prepare_sector_network import co2_emissions_year
from helper import update_config_with_sector_opts

#consolidate and rename
def rename_techs(label):

    prefix_to_remove = [
        "residential ",
        "services ",
        "urban ",
        "rural ",
        "central ",
        "decentral "
    ]

    rename_if_contains = [
        "CHP",
        "gas boiler",
        "biogas",
        "solar thermal",
        "air heat pump",
        "ground heat pump",
        "resistive heater",
        "Fischer-Tropsch"
    ]

    rename_if_contains_dict = {
        "water tanks": "hot water storage",
        "retrofitting": "building retrofitting",
        # "H2 Electrolysis": "hydrogen storage",
        # "H2 Fuel Cell": "hydrogen storage",
        # "H2 pipeline": "hydrogen storage",
        "battery": "battery storage",
        # "CC": "CC"
    }

    rename = {
        "solar": "solar PV",
        "Sabatier": "methanation",
        "offwind": "offshore wind",
        "offwind-ac": "offshore wind (AC)",
        "offwind-dc": "offshore wind (DC)",
        "onwind": "onshore wind",
        "ror": "hydroelectricity",
        "hydro": "hydroelectricity",
        "PHS": "hydroelectricity",
        "co2 Store": "DAC",
        "co2 stored": "CO2 sequestration",
        "AC": "transmission lines",
        "DC": "transmission lines",
        "B2B": "transmission lines"
    }

    for ptr in prefix_to_remove:
        if label[:len(ptr)] == ptr:
            label = label[len(ptr):]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old,new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old,new in rename.items():
        if old == label:
            label = new
    return label


preferred_order = pd.Index([
    "transmission lines",
    "hydroelectricity",
    "hydro reservoir",
    "run of river",
    "pumped hydro storage",
    "solid biomass",
    "biogas",
    "onshore wind",
    "offshore wind",
    "offshore wind (AC)",
    "offshore wind (DC)",
    "solar PV",
    "solar thermal",
    "solar rooftop",
    "solar",
    "building retrofitting",
    "ground heat pump",
    "air heat pump",
    "heat pump",
    "resistive heater",
    "power-to-heat",
    "gas-to-power/heat",
    "CHP",
    "OCGT",
    "gas boiler",
    "gas",
    "natural gas",
    "helmeth",
    "methanation",
    "hydrogen storage",
    "power-to-gas",
    "power-to-liquid",
    "battery storage",
    "hot water storage",
    "CO2 sequestration"
])

def plot_costs():


    cost_df = pd.read_csv(
        snakemake.input.costs,
        index_col=list(range(3)),
        header=list(range(n_header))
    )

    df = cost_df.groupby(cost_df.index.get_level_values(2)).sum()

    #convert to billions
    df = df / 1e9

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.max(axis=1) < snakemake.config['plotting']['costs_threshold']]

    print("dropping")

    print(df.loc[to_drop])

    df = df.drop(to_drop)

    print(df.sum())

    new_index = preferred_order.intersection(df.index).append(df.index.difference(preferred_order))

    new_columns = df.sum().sort_values().index

    fig, ax = plt.subplots(figsize=(12,8))

    df.loc[new_index,new_columns].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[snakemake.config['plotting']['tech_colors'][i] for i in new_index]
    )

    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([0,snakemake.config['plotting']['costs_max']])

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

    ax.grid(axis='x')

    ax.legend(handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1,1], frameon=False)

    fig.savefig(snakemake.output.costs, bbox_inches='tight')


def plot_energy():

    energy_df = pd.read_csv(
        snakemake.input.energy,
        index_col=list(range(2)),
        header=list(range(n_header))
    )

    df = energy_df.groupby(energy_df.index.get_level_values(1)).sum()

    #convert MWh to TWh
    df = df / 1e6

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.abs().max(axis=1) < snakemake.config['plotting']['energy_threshold']]

    print("dropping")

    print(df.loc[to_drop])

    df = df.drop(to_drop)

    print(df.sum())

    print(df)

    new_index = preferred_order.intersection(df.index).append(df.index.difference(preferred_order))

    new_columns = df.columns.sort_values()

    fig, ax = plt.subplots(figsize=(12,8))

    print(df.loc[new_index, new_columns])

    df.loc[new_index, new_columns].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[snakemake.config['plotting']['tech_colors'][i] for i in new_index]
    )

    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([snakemake.config['plotting']['energy_min'], snakemake.config['plotting']['energy_max']])

    ax.set_ylabel("Energy [TWh/a]")

    ax.set_xlabel("")

    ax.grid(axis="x")

    ax.legend(handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False)

    fig.savefig(snakemake.output.energy, bbox_inches='tight')



def plot_balances():

    co2_carriers = ["co2", "co2 stored", "process emissions"]

    balances_df = pd.read_csv(
        snakemake.input.balances,
        index_col=list(range(3)),
        header=list(range(n_header))
    )

    balances = {i.replace(" ","_"): [i] for i in balances_df.index.levels[0]}
    balances["energy"] = [i for i in balances_df.index.levels[0] if i not in co2_carriers]

    for k, v in balances.items():

        df = balances_df.loc[v]
        df = df.groupby(df.index.get_level_values(2)).sum()

        #convert MWh to TWh
        df = df / 1e6

        #remove trailing link ports
        df.index = [i[:-1] if ((i != "co2") and (i[-1:] in ["0","1","2","3"])) else i for i in df.index]

        df = df.groupby(df.index.map(rename_techs)).sum()

        to_drop = df.index[df.abs().max(axis=1) < snakemake.config['plotting']['energy_threshold']/10]

        print("dropping")

        print(df.loc[to_drop])

        df = df.drop(to_drop)

        print(df.sum())

        if df.empty:
            continue

        new_index = preferred_order.intersection(df.index).append(df.index.difference(preferred_order))

        new_columns = df.columns.sort_values()

        fig, ax = plt.subplots(figsize=(12,8))

        df.loc[new_index,new_columns].T.plot(kind="bar",ax=ax,stacked=True,color=[snakemake.config['plotting']['tech_colors'][i] for i in new_index])


        handles,labels = ax.get_legend_handles_labels()

        handles.reverse()
        labels.reverse()

        if v[0] in co2_carriers:
            ax.set_ylabel("CO2 [MtCO2/a]")
        else:
            ax.set_ylabel("Energy [TWh/a]")

        ax.set_xlabel("")

        ax.grid(axis="x")

        ax.legend(handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False)


        fig.savefig(snakemake.output.balances[:-10] + k + ".pdf", bbox_inches='tight')


def historical_emissions(cts):
    """
    read historical emissions to add them to the carbon budget plot
    """
    #https://www.eea.europa.eu/data-and-maps/data/national-emissions-reported-to-the-unfccc-and-to-the-eu-greenhouse-gas-monitoring-mechanism-16
    #downloaded 201228 (modified by EEA last on 201221)
    fn = "data/eea/UNFCCC_v23.csv"
    df = pd.read_csv(fn, encoding="latin-1")
    df.loc[df["Year"] == "1985-1987","Year"] = 1986
    df["Year"] = df["Year"].astype(int)
    df = df.set_index(['Year', 'Sector_name', 'Country_code', 'Pollutant_name']).sort_index()

    e = pd.Series()
    e["electricity"] = '1.A.1.a - Public Electricity and Heat Production'
    e['residential non-elec'] = '1.A.4.b - Residential'
    e['services non-elec'] = '1.A.4.a - Commercial/Institutional'
    e['rail non-elec'] = "1.A.3.c - Railways"
    e["road non-elec"] = '1.A.3.b - Road Transportation'
    e["domestic navigation"] = "1.A.3.d - Domestic Navigation"
    e['international navigation'] = '1.D.1.b - International Navigation'
    e["domestic aviation"] = '1.A.3.a - Domestic Aviation'
    e["international aviation"] = '1.D.1.a - International Aviation'
    e['total energy'] = '1 - Energy'
    e['industrial processes'] = '2 - Industrial Processes and Product Use'
    e['agriculture'] = '3 - Agriculture'
    e['LULUCF'] = '4 - Land Use, Land-Use Change and Forestry'
    e['waste management'] = '5 - Waste management'
    e['other'] = '6 - Other Sector'
    e['indirect'] = 'ind_CO2 - Indirect CO2'
    e["total wL"] = "Total (with LULUCF)"
    e["total woL"] = "Total (without LULUCF)"

    pol = ["CO2"] # ["All greenhouse gases - (CO2 equivalent)"]
    cts
    if "GB" in cts:
        cts.remove("GB")
        cts.append("UK")

    year = np.arange(1990,2018).tolist()

    idx = pd.IndexSlice
    co2_totals = df.loc[idx[year,e.values,cts,pol],"emissions"].unstack("Year").rename(index=pd.Series(e.index,e.values))

    co2_totals = (1/1e6)*co2_totals.groupby(level=0, axis=0).sum() #Gton CO2

    co2_totals.loc['industrial non-elec'] = co2_totals.loc['total energy'] - co2_totals.loc[['electricity', 'services non-elec','residential non-elec', 'road non-elec',
                                                                              'rail non-elec', 'domestic aviation', 'international aviation', 'domestic navigation',
                                                                              'international navigation']].sum()

    emissions = co2_totals.loc["electricity"]
    if "T" in opts:
        emissions += co2_totals.loc[[i+ " non-elec" for i in ["rail","road"]]].sum()
    if "H" in opts:
        emissions += co2_totals.loc[[i+ " non-elec" for i in ["residential","services"]]].sum()
    if "I" in opts:
        emissions += co2_totals.loc[["industrial non-elec","industrial processes",
                                          "domestic aviation","international aviation",
                                          "domestic navigation","international navigation"]].sum()
    return emissions



def plot_carbon_budget_distribution(input_eurostat):
    """
    Plot historical carbon emissions in the EU and decarbonization path
    """

    import matplotlib.gridspec as gridspec
    import seaborn as sns; sns.set()
    sns.set_style('ticks')
    plt.style.use('seaborn-ticks')
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20

    plt.figure(figsize=(10, 7))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs1[0,0])
    ax1.set_ylabel('CO$_2$ emissions (Gt per year)',fontsize=22)
    ax1.set_ylim([0,5])
    ax1.set_xlim([1990,snakemake.config['scenario']['planning_horizons'][-1]+1])

    path_cb = snakemake.config['results_dir'] + snakemake.config['run'] + '/csvs/'
    countries = pd.read_csv(snakemake.input.country_codes, index_col=1)
    cts = countries.index.to_list()
    e_1990 = co2_emissions_year(cts, input_eurostat, opts, year=1990)
    CO2_CAP=pd.read_csv(path_cb + 'carbon_budget_distribution.csv',
                        index_col=0)


    ax1.plot(e_1990*CO2_CAP[o],linewidth=3,
                         color='dodgerblue', label=None)

    emissions = historical_emissions(cts)

    ax1.plot(emissions, color='black', linewidth=3, label=None)

    #plot commited and uder-discussion targets
    #(notice that historical emissions include all countries in the
    # network, but targets refer to EU)
    ax1.plot([2020],[0.8*emissions[1990]],
                     marker='*', markersize=12, markerfacecolor='black',
                     markeredgecolor='black')

    ax1.plot([2030],[0.45*emissions[1990]],
                     marker='*', markersize=12, markerfacecolor='white',
                     markeredgecolor='black')

    ax1.plot([2030],[0.6*emissions[1990]],
                     marker='*', markersize=12, markerfacecolor='black',
                     markeredgecolor='black')

    ax1.plot([2050, 2050],[x*emissions[1990] for x in [0.2, 0.05]],
                  color='gray', linewidth=2, marker='_', alpha=0.5)

    ax1.plot([2050],[0.01*emissions[1990]],
                     marker='*', markersize=12, markerfacecolor='white',
                     linewidth=0, markeredgecolor='black',
                     label='EU under-discussion target', zorder=10,
                     clip_on=False)

    ax1.plot([2050],[0.125*emissions[1990]],'ro',
                     marker='*', markersize=12, markerfacecolor='black',
                     markeredgecolor='black', label='EU commited target')

    ax1.legend(fancybox=True, fontsize=18, loc=(0.01,0.01),
                       facecolor='white', frameon=True)

    path_cb_plot = snakemake.config['results_dir'] + snakemake.config['run'] + '/graphs/'
    plt.savefig(path_cb_plot+'carbon_budget_plot.pdf', dpi=300)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('plot_summary')


    n_header = 4

    plot_costs()

    plot_energy()

    plot_balances()

    for sector_opts in snakemake.config['scenario']['sector_opts']:
        opts=sector_opts.split('-')
        for o in opts:
            if "cb" in o:
                plot_carbon_budget_distribution(snakemake.input.eurostat)
