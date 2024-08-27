# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates plots from summary CSV files.
"""

import logging

import matplotlib.pyplot as plt
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from plot_summary import rename_techs, preferred_order

logger = logging.getLogger(__name__)
plt.style.use("ggplot")
plt.rcParams['axes.facecolor'] = 'white'  # Background of the plot
plt.rcParams['figure.facecolor'] = 'white'  # Background of the figure

plt.rcParams['grid.color'] = 'gray'  # Gridline color
plt.rcParams['grid.linestyle'] = '-'  # Solid gridline
plt.rcParams['grid.alpha'] = 0.7  # Gridline transparency
plt.rcParams['grid.linewidth'] = 0.8  # Gridline width


grouper = {
    "CCGT": "power",
    "OCGT": "power",
    "lignite": "power",
    "coal": "power",
    'oil': "power",
    'nuclear': "power",
    "agriculture machinery oil": "agriculture",
    'coal for industry': "industry",
    'gas for industry': "industry",
    'gas for industry CC': "industry CC",
    'industry methanol': "industry",
    'process emissions': "industry",
    'naphtha for industry': "industry",
    'process emissions CC': "industry CC",
    "solid biomass for industry CC": "industry CC",
    "co2": "co2 atmosphere",
    'kerosene for aviation': "aviation",
    'shipping methanol': "shipping",
    'shipping oil': "shipping",
    'urban central gas CHP': "heat",
    'urban central gas boiler': "heat",
    'urban decentral gas boiler': "heat",
    'urban central gas CHP CC': "heat CC",
    'urban central solid biomass CHP CC': "heat CC",
    'urban decentral oil boiler': "heat",
    'rural gas boiler': "heat",
    'rural oil boiler': "heat",
    "land transport oil light": "land transport",
    "land transport oil heavy": "land transport",
    }

def plot_costs(cost_df, drop=None):


    df = cost_df.groupby(cost_df.index.get_level_values(2)).sum()
    
    
    # plot transport
    df_filtered = df[df.index.str.contains("land transport")].droplevel([1, 2, 3], axis=1)
    
    # Reshape the dataframe using stack to move the planning horizon into the rows
    df_stacked = df_filtered.stack(level=1)
    
    planning_horizons = df.columns.get_level_values('planning_horizon').unique().sort_values()
    scenarios = df_filtered.columns.get_level_values(0).unique()
    
    
    fig, axes = plt.subplots(
    nrows=1, ncols=len(planning_horizons), 
    figsize=(12, 8), 
    sharey=True  # This ensures that all subplots share the same y-axis
    )
    
    if len(planning_horizons) == 1:
        axes = [axes]
    
    for ax, year in zip(axes, planning_horizons):
        subset = df_filtered.xs(year, level='planning_horizon', axis=1) 
        # convert to million
        subset /= 1e6
        
        subset.T.plot(
            kind="bar",
            ax=ax,
            stacked=True,
            legend=False,
            color=[snakemake.config["plotting"]["tech_colors"][i] for i in subset.index]
        )
    
        # Set title and x-label
        ax.set_title(year)
        
        ax.set_xlabel("")
        
        # Set x-ticks as scenario names (level=0)
        ax.set_xticks(range(len(scenarios)))
        
        if ax == axes[0]:
            ax.set_ylabel("System Cost [EUR million per year]")
            
        
        ax.grid(axis="x")
        
        fig.savefig(snakemake.output.costs.split(".pdf")[0]+"-transport.pdf",
                    bbox_inches="tight")



    # convert to billions
    df = df / 1e9

    df = df.groupby(df.index.map(rename_techs)).sum()

    to_drop = df.index[df.max(axis=1) < snakemake.config["plotting"]["costs_threshold"]]

    logger.info(
        f"Dropping technology with costs below {snakemake.config['plotting']['costs_threshold']} EUR billion per year"
    )
    logger.debug(df.loc[to_drop])

    df = df.drop(to_drop)

    logger.info(f"Total system cost of {round(df.sum().iloc[0])} EUR billion per year")

    new_index = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )

    # new_columns = df.sum().sort_values().index
    if drop!=None:
        df = df.droplevel([1,2,3], axis=1)
        
    new_cols = df.T.sort_index(level="planning_horizon").index
    
    planning_horizons = df.columns.get_level_values('planning_horizon').unique().sort_values()
    scenarios = df.columns.get_level_values(0).unique()
    
    fig, axes = plt.subplots(
    nrows=1, ncols=len(planning_horizons), 
    figsize=(12, 8), 
    sharey=True  # This ensures that all subplots share the same y-axis
    )
    
    if len(planning_horizons) == 1:
        axes = [axes]
    
    for ax, year in zip(axes, planning_horizons):
        subset = df.xs(year, level='planning_horizon', axis=1) 
        
        subset.T[new_index].plot(
            kind="bar",
            ax=ax,
            stacked=True,
            legend=False,
            color=[snakemake.config["plotting"]["tech_colors"][i] for i in new_index]
        )
    
        # Set title and x-label
        ax.set_title(year)
        
        ax.set_xlabel("")
        
        # Set x-ticks as scenario names (level=0)
        ax.set_xticks(range(len(scenarios)))
        
        if ax == axes[0]:
            ax.set_ylabel("System Cost [EUR billion per year]")
            
        
        ax.grid(axis="x")
        
        ax.set_ylim([0, snakemake.config['plotting']["costs_max"]])
    

    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()
        
    axes[-1].legend(
        handles,
        labels,
        ncol=1,
        loc="upper left",
        bbox_to_anchor=[1, 1],
        frameon=False,
    )
    
    plt.tight_layout()
    
    fig.savefig(snakemake.output.costs,
                bbox_inches="tight")
    
    
    df.sum().unstack().T.plot()
    plt.ylabel("System Cost [EUR billion per year]")
    plt.xlabel("planning horizon")
    plt.legend(bbox_to_anchor=(1,1))
    plt.savefig(snakemake.output.costs.split(".pdf")[0]+"-total.pdf",
                bbox_inches="tight")
    


def plot_balances(balances_df, drop=None):
    
    co2_carriers = ["co2", "co2 stored", "process emissions"]
    balances = {i.replace(" ", "_"): [i] for i in balances_df.index.levels[0]}
    balances["energy"] = [
        i for i in balances_df.index.levels[0] if i not in co2_carriers
    ]

    for k, v in balances.items():
        df = balances_df.loc[v]
        df = df.groupby(df.index.get_level_values(2)).sum()

        # convert MWh to TWh
        df = df / 1e6

        # remove trailing link ports
        df.index = [
            (
                i[:-1]
                if (
                    (i not in ["co2", "NH3", "H2"])
                    and (i[-1:] in ["0", "1", "2", "3", "4"])
                )
                else i
            )
            for i in df.index
        ]
        
        if k =="co2":
            
            co2_b = df.rename(index=grouper).groupby(level=0).sum()
            co2_b = co2_b.droplevel([1,2,3], axis=1)
            co2_b.loc["co2 atmosphere"] = abs(co2_b.loc["co2 atmosphere"])*-1
            
            scenarios = df.columns.get_level_values(0).unique()
            
            fig, axes = plt.subplots(
                nrows=len(scenarios), ncols=1, 
            figsize=(8, 12), 
            sharey=True,
            sharex=True
            )
            
            
            for ax, scenario in zip(axes, scenarios):
                co2_b[scenario].T.plot(kind="area", stacked=True,
                                       title=scenario,
                                       ax=ax, legend=False)
                
                ax.set_xlabel("")
                ax.set_ylabel("")
            axes[2].set_ylabel("CO2 [MtCO2/a]")
            axes[2].legend(bbox_to_anchor=(1,1))
            
            fig.savefig(snakemake.output.balances[:-19] + "co2-stacked-area.pdf",
                        bbox_inches="tight")
            
                
            for carrier in co2_b.index:
                co2_b.loc[carrier].unstack().T.plot(title=carrier)
                plt.ylabel("CO2 [MtCO2/a]")
                plt.xlabel("")
                plt.savefig(snakemake.output.balances[:-19] + f"co2-{carrier}.pdf",
                            bbox_inches="tight")
        
            import seaborn as sns
            scenario1 = [scenarios[0]]
            diff = (co2_b.stack()-co2_b.stack()[scenario1].values)
            bool_index = abs(diff.groupby(level=0).sum().sum(axis=1))>20
            # Calculate the global min and max for the colormap
            global_min = diff.min().min()
            global_max = diff.max().max()
                
            fig, axes = plt.subplots(
                nrows=len(scenarios)-1, ncols=1, 
                figsize=(8, 12), sharex=True)
            i = 0
            for scenario in scenarios:
                if scenario in scenario1: continue
                        
                diff = co2_b[scenario] - co2_b[scenario1].values
                
                heatmap_data = diff.loc[bool_index]
    
                plt.figure(figsize=(12, 8))
                sns.heatmap(heatmap_data, cmap='coolwarm',
                            annot=True, fmt=".0f", linewidths=.5,
                            center=0, ax=axes[i],
                            vmin=global_min, vmax=global_max,
                            cbar_kws={'label': 'Difference in CO$_2$ emissions [MtCO$_2$/a]'}
                            )
                axes[i].set_title(f'Difference in Emissions ({scenario} vs {scenario1[0]})')
                axes[i].set_xlabel('')
                i += 1
            fig.savefig(snakemake.output.balances[:-19] + "co2-heatmap.pdf",
                        bbox_inches="tight")

            
     

        df = df.groupby(df.index.map(rename_techs)).sum()

        to_drop = df.index[
            df.abs().max(axis=1) < snakemake.config["plotting"]["energy_threshold"] / 10
        ]

        units = "MtCO2/a" if v[0] in co2_carriers else "TWh/a"
        
                        
        logger.debug(
            f"Dropping technology energy balance smaller than {snakemake.config['plotting']['energy_threshold']/10} {units}"
        )
        logger.debug(df.loc[to_drop])

        df = df.drop(to_drop)

        logger.debug(
            f"Total energy balance for {v} of {round(df.sum().iloc[0],2)} {units}"
        )

        if df.empty:
            continue
        
        
            
            
            
        new_index = preferred_order.intersection(df.index).append(
            df.index.difference(preferred_order)
        )
        
        if drop!=None:
            df = df.droplevel([1,2,3], axis=1)

        
        planning_horizons = df.columns.get_level_values('planning_horizon').unique().sort_values()
        scenarios = df.columns.get_level_values(0).unique()
        
        fig, axes = plt.subplots(
        nrows=1, ncols=len(planning_horizons), 
        figsize=(12, 8), 
        sharey=True  # This ensures that all subplots share the same y-axis
        )
        
        if len(planning_horizons) == 1:
            axes = [axes]
        
        for ax, year in zip(axes, planning_horizons):
            subset = df.xs(year, level='planning_horizon', axis=1) 
            
            subset.T[new_index].plot(
                kind="bar",
                ax=ax,
                stacked=True,
                legend=False,
                color=[snakemake.config["plotting"]["tech_colors"][i] for i in new_index]
            )
        
            # Set title and x-label
            ax.set_title(year)
            
            # Set x-ticks as scenario names (level=0)
            ax.set_xticks(range(len(scenarios)))
            
            if ax == axes[0]:
                if v[0] in co2_carriers:
                    ax.set_ylabel("CO2 [MtCO2/a]")
                else:
                    ax.set_ylabel("Energy [TWh/a]")
            
            ax.grid(axis="x")
        

        handles, labels = ax.get_legend_handles_labels()
        handles.reverse()
        labels.reverse()
        
        axes[-1].legend(
            handles,
            labels,
            ncol=1,
            loc="upper left",
            bbox_to_anchor=[1, 1],
            frameon=False,
        )
        
        plt.tight_layout()
        
        fig.savefig(snakemake.output.balances[:-10] + k + ".pdf", bbox_inches="tight")
   
        if k in [ 'land_transport_demand_heavy',  'land_transport_demand_light']:
            share = abs(df.div(df.loc[v].values, axis=1).drop(v)*100)
            
            years = share.columns.levels[1]
            fig, axes = plt.subplots(nrows=len(years), ncols=1,
                                     sharey=True, sharex=True,
                                     figsize=(5,8))
            for i, year in enumerate(years):
                share.xs(year, level=1, axis=1).T.plot(kind="bar",
                                                       stacked=True,
                                                       color=[snakemake.config["plotting"]["tech_colors"][i] for i in share.index],
                                                       title=year,
                                                       ax=axes[i],
                                                       legend=False)
                axes[i].set_ylabel("share [%]")
                for container in axes[i].containers:
                    labels = [f'{round(v.get_height())}%' if int(v.get_height()) != 0 else '' for v in container]
                    axes[i].bar_label(container, labels=labels,
                                      label_type='center')  # You can also use 'edge' or 'center'
                fig.savefig(snakemake.output.balances[:-19] + f"share-bar-{k}.pdf", bbox_inches="tight")
                
            
            line_styles = ["-", "--", ":", "-."]

            
            # If there are more scenarios than line styles, this will repeat the line styles
            line_styles = (line_styles * (len(share.columns.levels[0]) // len(line_styles) + 1))[:len(share.columns.levels[0])]

            
            fig, axes = plt.subplots(
                        nrows=1, ncols=len(share.index), 
                        figsize=(14, 5), 
                        sharey=True  # This ensures that all subplots share the same y-axis
                        )
            for i in range(len(share.index)):
                share.iloc[i].unstack().T.plot(title=share.index[i], ax=axes[i],
                                               legend=False,
                                               style=line_styles)
                axes[i].set_xlabel("")
                axes[i].grid(axis="x")
            axes[0].set_ylabel("share [%]")
            
            axes[1].legend(
                ncol=1,
                loc="upper left",
                frameon=False,
                )
            
            fig.savefig(snakemake.output.balances[:-19] + f"share-{k}.pdf", bbox_inches="tight")
            

def plot_prices(prices, drop=True):
    if drop:
        prices = prices.droplevel([1,2,3], axis=1)
    
    for carrier in prices.index: #["co2", "AC", "low voltage", "H2"]:
        fig, axes = plt.subplots()
        if carrier in ["co2", "process emissions"]:
            label = "CO$_2$ price [Eur/t]"
            sign = -1
        else:
            label = "price [Eur/MWh]"
            sign = 1
        (prices.loc[carrier].unstack().T*sign).plot(ax=axes, title=carrier)
        axes.set_ylabel(label)
        axes.set_xlabel("")
        axes.set_xlim(left=0)
        axes.grid(axis="x")
        axes.legend(
            ncol=1,
            # loc="lower right",
            frameon=False,
            )
        fig.savefig(snakemake.output.balances[:-19] + f"prices-{carrier}.pdf",
                    bbox_inches="tight")
#%%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("plot_summary_all")

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    
    n_header = 4
    
    prefix = snakemake.config["run"]["prefix"]
    path = snakemake.output[0].split("graphs")[0]
    scenarios = snakemake.config["run"]["name"]
    
    costs = {}
    balances = {}
    prices = {}
    for scenario in scenarios:
        try:
            costs[scenario] = pd.read_csv(f"{path}/{scenario}/csvs/costs.csv",
                                          index_col=list(range(3)),
                                          header=list(range(n_header)))
            balances[scenario] = pd.read_csv(f"{path}/{scenario}/csvs/supply_energy.csv",
                                          index_col=list(range(3)),
                                          header=list(range(n_header)))
            prices[scenario] = pd.read_csv(f"{path}/{scenario}/csvs/prices.csv",
                                                      index_col=[0],
                                                      header=list(range(n_header)))
        except FileNotFoundError:
            logger.info(f"{scenario} not solved yet.")
            
    costs = pd.concat(costs, axis=1)
    balances = pd.concat(balances, axis=1)
    prices = pd.concat(prices, axis=1)
    
    costs.to_csv(snakemake.output.costs_csv)
    balances.to_csv(snakemake.output.balances_csv)
    prices.to_csv(snakemake.output.prices_csv)
      
    plot_costs(costs, drop=True)

    plot_balances(balances, drop=True)
    
    plot_prices(prices, drop=True)


