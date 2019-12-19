


import pandas as pd

#allow plotting without Xwindows
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt



#consolidate and rename
def rename_techs(label):

    prefix_to_remove = ["residential ","services ","urban ","rural ","central ","decentral "]

    rename_if_contains = ["CHP","gas boiler","biogas","solar thermal","air heat pump","ground heat pump","resistive heater","Fischer-Tropsch"]

    rename_if_contains_dict = {"water tanks" : "hot water storage",
                               "retrofitting" : "building retrofitting",
                               "H2" : "hydrogen storage",
                               "battery" : "battery storage",
                               "CCS" : "CCS"}

    rename = {"solar" : "solar PV",
              "Sabatier" : "methanation",
              "offwind" : "offshore wind",
              "offwind-ac" : "offshore wind (AC)",
              "offwind-dc" : "offshore wind (DC)",
              "onwind" : "onshore wind",
              "ror" : "hydroelectricity",
              "hydro" : "hydroelectricity",
              "PHS" : "hydroelectricity",
              "co2 Store" : "DAC",
              "co2 stored" : "CO2 sequestration",
              "AC" : "transmission lines",
              "DC" : "transmission lines",
              "B2B" : "transmission lines"}

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


preferred_order = pd.Index(["transmission lines","hydroelectricity","hydro reservoir","run of river","pumped hydro storage","solid biomass","biogas","onshore wind","offshore wind","offshore wind (AC)","offshore wind (DC)","solar PV","solar thermal","solar","building retrofitting","ground heat pump","air heat pump","heat pump","resistive heater","power-to-heat","gas-to-power/heat","CHP","OCGT","gas boiler","gas","natural gas","helmeth","methanation","hydrogen storage","power-to-gas","power-to-liquid","battery storage","hot water storage","CO2 sequestration"])

def plot_costs():


    cost_df = pd.read_csv(snakemake.input.costs,index_col=list(range(3)),header=[0,1,2])


    df = cost_df.groupby(cost_df.index.get_level_values(2)).sum()

    #convert to billions
    df = df/1e9

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.max(axis=1) < snakemake.config['plotting']['costs_threshold']]

    print("dropping")

    print(df.loc[to_drop])

    df = df.drop(to_drop)

    print(df.sum())

    new_index = (preferred_order&df.index).append(df.index.difference(preferred_order))

    new_columns = df.sum().sort_values().index

    fig, ax = plt.subplots()
    fig.set_size_inches((12,8))

    df.loc[new_index,new_columns].T.plot(kind="bar",ax=ax,stacked=True,color=[snakemake.config['plotting']['tech_colors'][i] for i in new_index])


    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([0,snakemake.config['plotting']['costs_max']])

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

    ax.grid(axis="y")

    ax.legend(handles,labels,ncol=4,loc="upper left")


    fig.tight_layout()

    fig.savefig(snakemake.output.costs,transparent=True)


def plot_energy():

    energy_df = pd.read_csv(snakemake.input.energy,index_col=list(range(2)),header=[0,1,2])

    df = energy_df.groupby(energy_df.index.get_level_values(1)).sum()

    #convert MWh to TWh
    df = df/1e6

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.abs().max(axis=1) < snakemake.config['plotting']['energy_threshold']]

    print("dropping")

    print(df.loc[to_drop])

    df = df.drop(to_drop)

    print(df.sum())

    print(df)

    new_index = (preferred_order&df.index).append(df.index.difference(preferred_order))

    new_columns = df.columns.sort_values()

    fig, ax = plt.subplots()
    fig.set_size_inches((12,8))

    print(df.loc[new_index,new_columns])

    df.loc[new_index,new_columns].T.plot(kind="bar",ax=ax,stacked=True,color=[snakemake.config['plotting']['tech_colors'][i] for i in new_index])


    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([snakemake.config['plotting']['energy_min'],snakemake.config['plotting']['energy_max']])

    ax.set_ylabel("Energy [TWh/a]")

    ax.set_xlabel("")

    ax.grid(axis="y")

    ax.legend(handles,labels,ncol=4,loc="upper left")


    fig.tight_layout()

    fig.savefig(snakemake.output.energy,transparent=True)


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils import Dict
        import yaml
        snakemake = Dict()
        with open('config.yaml') as f:
            snakemake.config = yaml.load(f)
        snakemake.input = Dict()
        snakemake.output = Dict()

        for item in ["costs","energy"]:
            snakemake.input[item] = snakemake.config['summary_dir'] + '/{name}/csvs/{item}.csv'.format(name=snakemake.config['run'],item=item)
            snakemake.output[item] = snakemake.config['summary_dir'] + '/{name}/graphs/{item}.pdf'.format(name=snakemake.config['run'],item=item)

    plot_costs()

    plot_energy()
