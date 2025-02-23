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
    'offshore wind (AC)': "power",
    'offshore wind (DC)': "power",
    'offshore wind (Float)':"power",
    'onshore wind': "power",
    'solar PV': "power",
    'uranium': "power",
    'nuclear': "power",
    "agriculture machinery oil": "agriculture",
    'coal for industry': "industry",
    'gas for industry': "industry",
    'gas for industry CC': "industry CC",
    'industry methanol': "industry",
    'process emissions': "industry",
    'naphtha for industry': "industry",
    'solid biomass for industry': "industry",
    'process emissions CC': "industry CC",
    "solid biomass for industry CC": "industry CC",
    "co2": "co2 atmosphere",
    'kerosene for aviation': "aviation",
    'shipping methanol': "shipping",
    'shipping oil': "shipping",
    'urban central gas CHP': "heat",
    'urban central gas boiler': "heat",
    'urban decentral gas boiler': "heat",
    'resistive heater': "heat",
    'biomass boiler': "heat",
    'hot water storage': "heat",
    'air heat pump': "heat",
    'ground heat pump': "heat",
    'oil boiler': "heat",
    'gas boiler': "heat",
    'urban central gas CHP CC': "heat CC",
    'urban central solid biomass CHP CC': "heat CC",
    'urban decentral oil boiler': "heat",
    'rural gas boiler': "heat",
    'rural oil boiler': "heat",
    "land transport oil light": "land transport",
    "land transport oil heavy": "land transport",
    'BEV charger heavy': "land transport",
    'BEV charger light': "land transport",
    'CHP': 'heat',
    'H2': "H2",
    'H2 Electrolysis':"H2",
    'H2 Fuel Cell':"H2",
    'H2 pipeline':"H2",
    'SMR':"H2",
    'SMR CC':"H2",
    'battery storage': "power",
    'solar rooftop':"power",
    'solar thermal':"power",
    'solar-hsat':"power",
    'transmission lines': "power",
    'electricity distribution grid': "power",
    'hydroelectricity': "power",
    'land transport EV heavy': "land transport",
    'land transport EV light': "land transport",
    'land transport fuel cell': "land transport",
    'land transport oil': "land transport",
    'methanolisation': "methanol",
    'biogas': "biomass",
    'solid biomass': "biomass",
    'co2 sequestered':'CO2 sequestration',
    'methanation': "gas",
    
    }


def get_unchagend_fleet():
    import yaml 
    from prepare_sector_network import get
    
    transport_types = ["light", "heavy"]
    ref_year = 2024
    
    registrations = pd.read_csv("/home/lisa/Documents/playground/pypsa-eur/resources/test-after-merge-v2/car_registration_s_39.csv",
                                index_col=[0,1])
    historical_ev_share = pd.read_csv("/home/lisa/Documents/playground/pypsa-eur/resources/test-after-merge-v2/historical_ev_share_s_39.csv",
                                      index_col=[0], header=[0,1])
    
    shares = {}
    for scenario in scenarios:
        fn = path + scenario + "/configs/config.base_s_39_lv1.0___2025.yaml"
        with open(fn, 'r') as file:
            config = yaml.safe_load(file)
        years = config["scenario"]["planning_horizons"]
        for transport_type in transport_types:
            for year in years:
    
                factor = get(config["sector"]["car_reg_factor"][transport_type], year)
                
                reg = registrations.loc[transport_type].iloc[:,0] * factor
                today_ev = historical_ev_share.ffill().iloc[-1].loc[transport_type]
                
                changed_fleet = (1-(reg*(year-ref_year))-today_ev).clip(lower=0)
                            
                shares[scenario, transport_type, year] = changed_fleet
                
        
def plot_costs(cost_df, drop=None):

    wished_scenarios = ["base", "fast", "slow", "high-demand", "low-demand", "mandate"]
    df = cost_df.groupby(cost_df.index.get_level_values(2)).sum()
    
    selector = df.columns.get_level_values(0).isin(wished_scenarios)
    df = df.loc[:, selector]
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
        subset /= 1e9
        
        subset.T.plot(
            kind="bar",
            ax=ax,
            stacked=True,
            # legend=False,
            color=[snakemake.config["plotting"]["tech_colors"][i] for i in subset.index]
        )
    
        # Set title and x-label
        ax.set_title(year)
        
        ax.set_xlabel("")
        
        # Set x-ticks as scenario names (level=0)
        ax.set_xticks(range(len(scenarios)))
        
        if ax == axes[0]:
            ax.set_ylabel("System Cost [EUR billion per year]")
            
        
        ax.grid(axis="x")
        
        fig.savefig(snakemake.output.costs.split(".pdf")[0]+"-transport.pdf",
                    bbox_inches="tight")



    # convert to billions
    df = df / 1e9

    df = df.groupby(df.index.map(rename_techs)).sum()
    
    df_grouped = df.rename(index=grouper).groupby(level=0).sum()

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
    figsize=(12, 11), 
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
        
        ax.set_ylim([0, 1.1*df.sum().max()])
    

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
    
    # together = df.sum().unstack().T
    # wished_scenarios = ["base", "fast", "slow"]
    # fig, ax = plt.subplots()
    # for sc in wished_scenarios:
    #     together[sc].plot(ax=ax)
    #     area = together.loc[:,together.columns.str.contains(sc)&~together.columns.str.contains("lock")
    #                         &~together.columns.str.contains("2")]
        
    #     lower_bound = area.min(axis=1)
    #     upper_bound = area.max(axis=1)
    #     ax.fill_between(area.index, lower_bound, upper_bound, alpha=0.2, label=f"{sc} range")
        
    # ax.set_ylabel("System Cost [EUR billion per year]")
    # ax.set_xlabel("planning horizon")
    


def split_positive_negative(df):
    """
    Splits each row in a DataFrame into two rows: one containing only positive values and the other only negative values.
    Both rows retain the same index as the original row.
    
    Args:
    df (pd.DataFrame): The input DataFrame with numeric values.
    
    Returns:
    pd.DataFrame: A DataFrame with rows split into positive and negative values.
    """
    positive_df = df.clip(lower=0)
    negative_df = df.clip(upper=0)
    
    # Append negative_df after positive_df, retaining the same index
    result_df = pd.concat([positive_df, negative_df])
    
    return result_df.sort_index()

def plot_balances(balances, drop=None):
    
    co2_carriers = ["co2", "co2 stored", "process emissions"]
    balances_df = {i.replace(" ", "_"): [i] for i in balances.index.levels[0]}
    balances_df["energy"] = [
        i for i in balances.index.levels[0] if i not in co2_carriers
    ]
    
    wished_scenarios = ["base", "fast", "slow", "high-demand", "low-demand"]

    for k, v in balances_df.items():
        df = balances.loc[v]
        df = df.groupby(df.index.get_level_values(2)).sum()
        
        # select only some scenarios
        selector = df.columns.get_level_values(0).isin(wished_scenarios)
        df = df.loc[:, selector]
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
        
        # if k=="AC":
        #     to_plot = df.droplevel([1,2,3], axis=1)[["base", "dsm", "v2g"]]
        #     diff = to_plot.subtract(to_plot.base, axis=0)
        #     gens = ["offwind-ac", "offwind-dc", "offwind-float", "onwind", "solar", "solar-hsat"]
            
        #     fig, ax = plt.subplots(nrows=1, ncols=len(gens), sharey=True,
        #                            figsize=(10,2.5))
        #     for col, gen in enumerate(gens):
        #         legend=True if col==0 else False
        #         diff.loc[gen].unstack().T.iloc[1:,1:].plot(kind="bar", title=gen,
        #                                                    ax=ax[col], legend=legend)
        #         ax[col].set_xlabel("")
        #     ax[0].set_ylabel("change in generation [TWh]")
        #     ax[2].set_xlabel("planning horizon")
            
        #     fig.savefig(snakemake.output.balances[:-19] + "gen-diff-dsm-v2g.pdf",
        #                 bbox_inches="tight")
        
            # dist = to_plot.loc["electricity distribution grid"].sum().unstack().T
            
            # caps = capacities.xs("electricity distribution grid", level=1).droplevel([1,2,3], axis=1)[["base", "dsm", "v2g"]].iloc[0,:].unstack().T
            
            # fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(8,3), sharex=True)
            # ((dist.divide(dist.base, axis=0)[["dsm", "v2g"]]-1)*100).plot(kind="bar", ax=ax[0], title="Difference in flow", legend=False)
            # ((caps.divide(caps.base, axis=0)[["dsm", "v2g"]]-1)*100).plot(kind="bar", ax=ax[1], title="Difference in installed capacity")
            # ax[1].set_xlabel("planning horizon")
            # ax[0].set_ylabel("Difference \n [%]")
            # ax[1].set_ylabel("Difference \n [%]")
            # fig.savefig(snakemake.output.balances[:-19] + "elec-dist-diff-dsm-v2g.pdf",
            #                 bbox_inches="tight")
            
            # # Compute the data to plot
            # to_plot = costs.sum().droplevel([1, 2, 3]).loc[["base", "dsm", "v2g"]].unstack().T / 1e9
            # percentage = ((to_plot.divide(to_plot.base, axis=0) - 1) * 100)
            
            # # Plot the bar chart
            # ax = to_plot.subtract(to_plot.base, axis=0).iloc[1:, 1:].plot(kind="bar")
            # plt.ylabel("Cost difference \n [bn Euro/a]")
            # plt.xlabel("Planning horizon")
            
            # # Add percentage text to each bar
            # j = 0
            # for container in ax.containers:  # Iterate through the bar groups
            #     for i, bar in enumerate(container):  # Iterate through individual bars
            #         x = bar.get_x() + bar.get_width() / 2  # X-coordinate (center of the bar)
            #         y = bar.get_height()  # Height of the bar
            #         if abs(y) > 1:  # Skip bars with zero height
            #             # Get the corresponding percentage value
            #             percentage_value = percentage.iloc[1:, 1:].iloc[i, j]
            #             # Annotate the bar
            #             ax.annotate(
            #                 f"{percentage_value:.1f}%",  # Format as percentage
            #                 (x, y),  # Position above the bar
            #                 ha="center", va="bottom", fontsize=9  # Styling
            #             )
            #         if i==3: j+=1

            
            # plt.ylabel("cost difference \n [bn Euro/a]")
            # plt.xlabel("planning horizon")
            # plt.savefig(snakemake.output.balances[:-19] + "cost-diff-dsm-v2g.pdf",
            #                 bbox_inches="tight")
            
        if k =="co2":
            
            co2_b = df.rename(index=grouper).groupby(level=0).sum()
            co2_b = co2_b.droplevel([1,2,3], axis=1)
            co2_b = split_positive_negative(co2_b)
            co2_b = co2_b[abs(co2_b.sum(axis=1))>5]
            #co2_b.loc["co2 atmosphere"] = abs(co2_b.loc["co2 atmosphere"])*-1
            
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
            
            grouped = co2_b.groupby(level=0).sum()
            for carrier in grouped.index:
                grouped.loc[carrier].unstack().T.plot(title=carrier)
                plt.ylabel("CO2 [MtCO2/a]")
                plt.xlabel("")
                plt.savefig(snakemake.output.balances[:-19] + f"co2-{carrier}.pdf",
                            bbox_inches="tight")
        
            import seaborn as sns
            scenarios = ["base", "fast", "slow"]
            if all([sc in co2_b.stack().columns for sc in scenarios]):
                scenario1 = [scenarios[0]]
                diff = (co2_b.stack()[scenarios]-co2_b.stack()[scenario1].values)
                diff = diff.groupby(level=[0,1]).sum()
                bool_index = (abs(diff)>2).any(axis=1).groupby(level=0).any()
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
                    diff = diff.groupby(level=0).sum()
                    heatmap_data = diff.loc[bool_index].drop("2025", axis=1)
        
                    plt.figure(figsize=(12, 8))
                    sns.heatmap(heatmap_data, cmap='coolwarm',
                                annot=True, fmt=".0f", linewidths=.5,
                                center=0, ax=axes[i],
                                vmin=global_min, vmax=global_max,
                                cbar=False,
                                #cbar_kws={'label': 'Difference in CO$_2$ emissions [MtCO$_2$/a]'}
                                )
                    axes[i].set_title(f'{scenario} vs {scenario1[0]}')
                    axes[i].set_xlabel('')
                    i += 1
                
                cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])  # Position of the colorbar [left, bottom, width, height]
                norm = plt.Normalize(vmin=global_min, vmax=global_max)
                sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=norm)
                sm.set_array([])  # Only needed for matplotlib < 3.1
                cbar = fig.colorbar(sm, cbar_ax, orientation='vertical')
                cbar.set_label('Difference in CO$_2$ emissions [MtCO$_2$/a]')
    
                fig.savefig(snakemake.output.balances[:-19] + "co2-heatmap.pdf",
                            bbox_inches="tight")
            
            scenarios = ["base", "low-demand", "high-demand"]
            if all([sc in co2_b.stack().columns for sc in scenarios]):
                scenario1 = [scenarios[0]]
                diff = (co2_b.stack()[scenarios]-co2_b.stack()[scenario1].values)
                diff = diff.groupby(level=[0,1]).sum()
                bool_index = (abs(diff)>2).any(axis=1).groupby(level=0).any()
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
                    diff = diff.groupby(level=0).sum()
                    heatmap_data = diff.loc[bool_index].drop("2025", axis=1)
        
                    plt.figure(figsize=(12, 8))
                    sns.heatmap(heatmap_data, cmap='coolwarm',
                                annot=True, fmt=".0f", linewidths=.5,
                                center=0, ax=axes[i],
                                vmin=global_min, vmax=global_max,
                                cbar=False,
                                #cbar_kws={'label': 'Difference in CO$_2$ emissions [MtCO$_2$/a]'}
                                )
                    axes[i].set_title(f'{scenario} vs {scenario1[0]}')
                    axes[i].set_xlabel('')
                    i += 1
                
                cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])  # Position of the colorbar [left, bottom, width, height]
                norm = plt.Normalize(vmin=global_min, vmax=global_max)
                sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=norm)
                sm.set_array([])  # Only needed for matplotlib < 3.1
                cbar = fig.colorbar(sm, cbar_ax, orientation='vertical')
                cbar.set_label('Difference in CO$_2$ emissions [MtCO$_2$/a]')
    
                fig.savefig(snakemake.output.balances[:-19] + "co2-heatmap-demand.pdf",
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
            
            exogen = shares.droplevel([1,2,3], axis=1)
            
            exogen_grouped = exogen[exogen.index.str.contains(k.split("_")[-1])].sum().unstack()
            
            exogen_share = exogen_grouped["existing"].div(exogen_grouped["demand"]).unstack().T
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
                if share.index[i] == "land transport oil":
                    a=exogen_share.rename(columns = lambda x: x + " exogen")
                    (a*100).plot(ax=axes[i], legend=False, style= "-.", color="black",
                                 linewidth=2)
                axes[i].set_xlabel("")
                axes[i].grid(axis="x")
            axes[0].set_ylabel("share [%]")
            
            axes[1].legend(
                ncol=1,
                loc="upper left",
                frameon=False,
                )
            
            fig.savefig(snakemake.output.balances[:-19] + f"share-{k}.pdf", bbox_inches="tight")

def plot_shares_area(balances, shares):

    rename_dict = {"EV": "BEV",
                   "fuel cell": "FCEV",
                   "oil": "ICE"}
    
    fig, axes = plt.subplots(
                nrows=2,
                ncols=3, 
                sharex=True,
                figsize=(14, 10), 
                sharey=True  # This ensures that all subplots share the same y-axis
                )
    carriers = ["land transport demand light"], ["land transport demand heavy"]
    for row, v in enumerate(carriers):
        df = balances.loc[v]
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
        
        df = df.groupby(df.index.map(rename_techs)).sum()

        to_drop = df.index[
            df.abs().max(axis=1) < snakemake.config["plotting"]["energy_threshold"] / 10
        ]

        units ="TWh/a"
        
                        
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
        
        
        df = df.droplevel([1,2,3], axis=1)

        
        planning_horizons = df.columns.get_level_values('planning_horizon').unique().sort_values()
        scenarios = df.columns.get_level_values(0).unique()
        
        share = abs(df.div(df.loc[v].values, axis=1).drop(v)*100)
        
        exogen = shares.droplevel([1,2,3], axis=1)
        exogen_grouped = exogen[exogen.index.str.contains(v[0].split("demand ")[-1])].sum().unstack()
        exogen_share = exogen_grouped["existing"].div(exogen_grouped["demand"]).unstack().T
        
           
        wished_scenarios = ["fast", "base", "slow"]
        color = ['g','b', 'r']
        line_styles = ["-", "--", ":", "-."]
        

        for i in range(len(share.index)):
            title = share.index[i].replace("land transport ", "").replace(" light", "") if row==0 else ""
            if title!="":
                title = rename_dict[title]
            share.iloc[i].reindex(wished_scenarios, level=0).unstack().T.plot(title=title,
                                                                              ax=axes[row, i],
                                           legend=False,
                                           style=line_styles[:len(wished_scenarios)],
                                           color=color)
            
            # Define color for the scenario's contour areas
            # color = axes[i].get_lines()[-1].get_color()  # Color of the last plotted line for consistency
            

            if share.index[i] == "land transport oil":
                a=(exogen_share.reindex(columns=wished_scenarios)
                   .dropna(how="all", axis=1)
                   .rename(columns = lambda x: x + " exogen"))
                (a*100).plot(ax=axes[row, i], legend=False, style=line_styles[:len(wished_scenarios)],
                             color="gray",
                             linewidth=2, alpha=0.5)
                
                axes[row, i].text(
                        x=a.index[-1],  # Position near the end of the x-axis
                        y=90, # (a * 100).iloc[-1].values[0],  # Position based on the last data point value
                        s="exogenous", 
                        color="gray", 
                        fontsize=10, 
                        ha="right",
                        va="center"
                    )
                
            for j, sc in enumerate(wished_scenarios):
               
                # smaller variation
                filter_b = (share.columns.get_level_values(0).str.contains(sc) & 
                ~share.columns.get_level_values(0).str.contains("2"))
                area = share.loc[:,filter_b].iloc[i]
                axes[row, i].fill_between(
                        area.unstack().columns,
                        area.unstack().min(),
                        area.unstack().max(),
                        color=color[j],
                        alpha=0.1,
                        label="+/-10% CAPEX" if (j==0 and row==0) else "",
                    )
                
                # larger variation
                area = share.loc[:,share.columns.get_level_values(0).str.contains(sc)].iloc[i]
                axes[row, i].fill_between(
                        area.unstack().columns,
                        area.unstack().min(),
                        area.unstack().max(),
                        color=color[j],
                        alpha=0.05,
                        label="+/-20% CAPEX" if (j==0 and row==0) else "",
                    )
                
            axes[row, i].set_xlabel("")
            axes[row, i].grid(axis="x")
        axes[row, 0].set_ylabel("share [%]")
        
        axes[0, 1].legend(
            ncol=1,
            loc="upper left",
            frameon=False,
            )

    # Add row titles for "light" and "heavy"
    fig.text(0.5, 0.95, "Light", ha='center', va='center', fontsize=16, weight='bold')
    fig.text(0.5, 0.5, "Heavy", ha='center', va='center', fontsize=16, weight='bold')
    fig.savefig(snakemake.output.balances[:-19] + "area-share.pdf",
                bbox_inches="tight")
    
    
    
    # new
    # Update rcParams for global font adjustments
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
    })

    fig, axes = plt.subplots(
                nrows=2,
                ncols=3, 
                sharex=True,
                figsize=(10, 5), 
                sharey=True  # This ensures that all subplots share the same y-axis
                )
    carriers = ["land transport demand light"], ["land transport demand heavy"]
    for row, v in enumerate(carriers):
        df = balances.loc[v]
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
        
        df = df.groupby(df.index.map(rename_techs)).sum()

        to_drop = df.index[
            df.abs().max(axis=1) < snakemake.config["plotting"]["energy_threshold"] / 10
        ]

        units ="TWh/a"
        
                        
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
        
        
        df = df.droplevel([1,2,3], axis=1)

        
        planning_horizons = df.columns.get_level_values('planning_horizon').unique().sort_values()
        scenarios = df.columns.get_level_values(0).unique()
        
        share = abs(df.div(df.loc[v].values, axis=1).drop(v)*100)
        
        exogen = shares.droplevel([1,2,3], axis=1)
        exogen_grouped = exogen[exogen.index.str.contains(v[0].split("demand ")[-1])].sum().unstack()
        exogen_share = exogen_grouped["existing"].div(exogen_grouped["demand"]).unstack().T
        
           
        wished_scenarios = ["fast", "base", "slow"]
        color = ['g','b', 'r']
        line_styles = ["-", "--", ":", "-."]
        

        for i in range(len(share.index)):
            title = share.index[i].replace("land transport ", "").replace(" light", "") if row==0 else ""
            if title!="":
                title = rename_dict[title]
            share.iloc[i].reindex(wished_scenarios, level=0).unstack().T.plot(title=title,
                                                                              ax=axes[row, i],
                                           legend=False,
                                           style=line_styles[:len(wished_scenarios)],
                                           color=color)
            
            # Define color for the scenario's contour areas
            # color = axes[i].get_lines()[-1].get_color()  # Color of the last plotted line for consistency
            

            if share.index[i] == "land transport oil":
                a=(exogen_share.reindex(columns=wished_scenarios)
                   .dropna(how="all", axis=1)
                   .rename(columns = lambda x: x + " exogen"))
                (a*100).plot(ax=axes[row, i], legend=False, style=line_styles[:len(wished_scenarios)],
                             color="gray",
                             linewidth=2, alpha=0.5)
                
                axes[row, i].text(
                        x=a.index[-1],  # Position near the end of the x-axis
                        y=90, # (a * 100).iloc[-1].values[0],  # Position based on the last data point value
                        s="exogenous", 
                        color="gray", 
                        fontsize=12, 
                        ha="right",
                        va="center"
                    )
                
            for j, sc in enumerate(wished_scenarios):
               
                # smaller variation
                # filter_b = (share.columns.get_level_values(0).str.contains(sc) & 
                # ~share.columns.get_level_values(0).str.contains("2"))
                # area = share.loc[:,filter_b].iloc[i]
                # axes[row, i].fill_between(
                #         area.unstack().columns,
                #         area.unstack().min(),
                #         area.unstack().max(),
                #         color=color[j],
                #         alpha=0.1,
                #         # label="+/-10% CAPEX" if (j==0 and row==0) else "",
                #     )
                
                # larger variation
                area = share.loc[:,share.columns.get_level_values(0).str.contains(sc)].iloc[i]
                axes[row, i].fill_between(
                        area.unstack().columns,
                        area.unstack().min(),
                        area.unstack().max(),
                        color=color[j],
                        alpha=0.05,
                        # label="+/-20% CAPEX" if (j==0 and row==0) else "",
                    )
                
            axes[row, i].set_xlabel("")
            axes[row, i].grid(axis="x")
        axes[row, 0].set_ylabel("share [%]")
        
        axes[0, 1].legend(
            ncol=1,
            loc="upper left",
            frameon=False,
            )

    # Add row titles for "light" and "heavy"
    fig.text(0.5, 0.99, "Light-duty", ha='center', va='center', fontsize=14, weight='bold')
    fig.text(0.5, 0.5, "Heavy-duty", ha='center', va='center', fontsize=14, weight='bold')
    fig.savefig(snakemake.output.balances[:-19] + "area-share-largevariations.pdf",
                bbox_inches="tight")
            

def plot_tech_assumptions():
    fn = "/home/lisa/Documents/own_projects/endogenous_transport/data/total_costs.csv"
    new_costs = pd.read_csv(fn, index_col=[0,1])
    
    car_keys = {'FCEV heavy': ['FCV Bus city',
      'FCV Coach',
      'FCV Truck Semi-Trailer max 50 tons',
      'FCV Truck Solo max 26 tons',
      'FCV Truck Trailer max 56 tons'],
     'ICE heavy': ['Diesel Bus city',
      'Diesel Coach',
      'Diesel Truck Semi-Trailer max 50 tons',
      'Diesel Truck Solo max 26 tons',
      'Diesel Truck Trailer max 56 tons'],
     'BEV heavy': ['BEV Bus city',
      'BEV Coach',
      'BEV Truck Semi-Trailer max 50 tons',
      'BEV Truck Solo max 26 tons',
      'BEV Truck Trailer max 56 tons'],
     'BEV light': ["Battery electric (passenger cars)"],
     'FCEV light': ["Hydrogen fuel cell (passenger cars)"],
     'ICE light': ["Liquid fuels ICE (passenger cars)"]
     }
    
    # Reverse the car_keys dictionary
    reverse_car_keys = {v: k for k, values in car_keys.items() for v in values}
    
    a = new_costs.stack().swaplevel(0,1)
    
    a = a.rename(index=reverse_car_keys, level=0).groupby(level=[0,1,2]).mean()

    a = a.loc[car_keys.keys()]
    
    to_keep = ['FOM', 'investment', 'fixed', 'efficiency', 'lifetime']
    
    a = a[a.index.get_level_values(2).isin(to_keep)]
    
    units = {"FOM": "% of investment/year",
             "investment": "EUR/vehicle",
             "fixed": "EUR/vehicle/a",
             "efficiency": "kWh/km",
             "lifetime": "years"}
    total = pd.DataFrame()
    total["value"] = a
    total["unit"] = a.reset_index()["level_2"].map(units).values
    
    # efficiencies light duty
    eff_light = {'BEV light': 0.188,
                 'FCEV light': 0.33, 
                 "ICE light": 0.622}
    # Filter rows where level=2 of the index is "efficiency"
    efficiency_rows = total.index.get_level_values(2) == "efficiency"
    
    # Iterate through eff_light and update matching rows in the DataFrame
    for tech, new_value in eff_light.items():
        # Replace value in total['value'] where technology matches and level=2 is "efficiency"
        mask = efficiency_rows & (total.index.get_level_values(0) == tech)
        total.loc[mask, 'value'] = new_value
    
def plot_prices(prices, drop=True):
    wished_scenarios = ["base", "fast", "slow", "high-demand", "low-demand", "mandate"]
    prices = prices[wished_scenarios]
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
        
def plot_capacities(capacities, drop=True):
    if drop:
        capacities = capacities.droplevel([1,2,3], axis=1)
    carriers = {"electrolysis": ["H2 Electrolysis"],
                "DAC": ["DAC"],
                "renewables": ['offwind-ac', 'offwind-dc', "offwind-float", 
                               "onwind", "solar", "solar rooftop", "solar-hsat"],
                "electric heating": ['urban decentral air heat pump',
                                     'urban central resistive heater',
                                     'urban central air heat pump',
                                     'urban decentral resistive heater',
                                     'rural ground heat pump',
                                     'rural air heat pump']
                }
    for carrier in carriers:
        df = capacities.droplevel(0).loc[carriers[carrier]]
        
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
            subset = df.xs(year, level='planning_horizon', axis=1)/1e3 
            subset = subset.rename(index=rename_techs).groupby(level=0).sum()
            
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
                ax.set_ylabel("Installed capacity [GW]")
                
            
            ax.grid(axis="x")

def plot_comparison(balances, capacities, drop=True):
    LHV_H2 = 33.33 # kWh/kg
    renewables = ['offwind-ac', 'offwind-dc', "offwind-float", 
                   "onwind", "solar", "solar rooftop", "solar-hsat"]
    wished_scenarios = ["fast", "slow", "low-demand", "high-demand", "mandate"]
    res = capacities.droplevel(0).loc[renewables].droplevel([1,2,3], axis=1)
    electrolysis = capacities.droplevel(0).loc["H2 Electrolysis"].droplevel([1,2,3])
    dac = balances.loc["co2"].droplevel(0).loc["DAC2"].droplevel([1,2,3])
    cost_df = costs.droplevel([0,1]).droplevel([1,2,3], axis=1).rename(index=rename_techs).groupby(level=0).sum()
    wished_scenarios = [sc for sc in wished_scenarios if sc in cost_df.columns.levels[0]]
    
    heat_balance = balances.loc[["rural heat", "urban decentral heat"]].droplevel([1]).groupby(level=1).sum().droplevel([1,2,3], axis=1).rename(index=rename_techs).groupby(level=0).sum()
    
    planning_horizons = res.columns.get_level_values('planning_horizon').unique().sort_values()
    
    nrows = 4
    ncols = len(planning_horizons[1:])
 
    fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols, 
    figsize=(12, 8), 
    sharex=True  
    )
    
    plt.subplots_adjust(
        left=0.1,  # Left margin
        right=0.9,  # Right margin
        top=0.9,  # Top margin
        bottom=0.1,  # Bottom margin
        hspace=0.7,  # Height (vertical) spacing between subplots
        wspace=0.3  # Width (horizontal) spacing between subplots
    )
    
    if len(planning_horizons[1:]) == 1:
        axes = [axes]
        
    # Manually set the y-axis sharing row-wise
    for i in range(nrows):
        for j in range(1, ncols):
            axes[i, j].sharey(axes[i, 0])  # Share y-axis of each column with the first column in each row
    
        
    for i, year in  enumerate(planning_horizons[1:]):
        
        # res capacities ------------------------------------------------------
        caps = (res - res["base"])[wished_scenarios].xs(year, level=1, axis=1).T/1e3
        for tech in ["offwind", "solar"]:
            offwind_i = caps.columns[caps.columns.str.contains(tech)]
            offshore = caps.loc[:,caps.columns.str.contains(tech)].sum(axis=1)
            caps.drop(offwind_i, axis=1, inplace=True)
            caps[tech]= offshore
        caps.rename(columns=rename_techs, inplace=True)
        # if i==0:
        #     legend=True
        # else:
        #     legend=False
        legend=False
        bars = caps.plot(kind="bar", stacked=True,
                  color = [snakemake.config["plotting"]["tech_colors"][i] for i in caps.T.index],
                  legend=legend, title=year, ax=axes[0, i])
        
        percent_change = res.sum().div(res["base"].sum())[wished_scenarios].xs(year, level=1)
        percent_change_labels = [f"{(p-1)*100:.0f}%" for p in percent_change]
        
        # Add text annotations on each bar
        j = 0
        for bar, label in zip(bars.patches, percent_change_labels):
            height = caps.sum(axis=1)[j]
            # Position the text slightly above the top of the bar
            axes[0, i].text(
                bar.get_x() + bar.get_width() / 2, 
                height + (0.02 * height),  # Adjust the vertical position
                label,
                ha='center', va='bottom',
                fontsize=10, color='black'
            )
            j += 1
            
        # electrified heating -----------------------------------------------
        # EU target is 13 million heat pumps sold per year by 2027
        typical_generation = 15 # MWh/a
        heat_techs = ["ground heat pump", "air heat pump"]#, "resistive heater"]
        heat_df = heat_balance.loc[heat_techs]/1e6
        number_pumps = (heat_df/typical_generation).xs(year, axis=1, level=1)
        diff_number = number_pumps[wished_scenarios].sub(number_pumps["base"], axis=0)
        to_plot = (heat_df - heat_df["base"]).xs(year, axis=1, level=1)[wished_scenarios]
        # to_plot.T.plot(kind="bar", stacked=True, ax=axes[1,i],
        #              color=[snakemake.config["plotting"]["tech_colors"][tech] for tech in to_plot.index],
        #              legend=legend)
        bars = diff_number.T.plot(kind="bar", stacked=True, ax=axes[1,i],
                          color=[snakemake.config["plotting"]["tech_colors"][tech] for tech in number_pumps.index],
                          legend=legend
                          )
        
        percent_change_labels = [f"{p:.0f} TWh" for p in to_plot.sum()]
        
        # Add text annotations on each bar
        j = 0
        for bar, label in zip(bars.patches, percent_change_labels):
            height = diff_number.sum()[j] #bar.get_height()
            # Position the text slightly above the top of the bar
            axes[1, i].text(
                bar.get_x() + bar.get_width() / 2, 
                height + (0.02 * height),  # Adjust the vertical position
                label,
                ha='center', va='bottom',
                fontsize=10, color='black'
            )
            j += 1
        

        caps = electrolysis.loc[wished_scenarios].xs(year, level=1)/1e3
        produced = balances.loc["H2"].droplevel(0).droplevel([1,2,3], axis=1).loc["H2 Electrolysis1"]
        # convert to million tonnes
        produced /= LHV_H2*1e6
        produced = produced.unstack().T
        diff = produced.subtract(produced.base, axis=0)
        bars = diff.loc[year, wished_scenarios].plot(kind="bar", legend=False,
                                        ax=axes[2,i],
                                        color=snakemake.config["plotting"]["tech_colors"]["H2"])
        if year == "2030":
            y_value =  20 - produced.loc["2030", "base"] 
            axes[2,i].axhline(y=y_value, ls="--",
                              color="black")
            # Add the text label near the line
            axes[2, i].text(
                x=-0.2,  # Position the text at the center of the x-axis (adjust as needed)
                y=y_value,  # Align the text vertically with the line
                s="RePowerEU",  # The text to display
                color="black",  # Text color
                ha='left',  # Horizontal alignment
                va='bottom',  # Vertical alignment, place text above the line
                fontsize=10  # Adjust the font size as needed
            )
            
        #     caps_labels = [f"{c:.0f} GW" for c in caps]
            
        #     # Add text annotations on each bar
        #     for bar, label in zip(bars.patches, caps_labels):
        #         height = bar.get_height()
        #         # Position the text slightly above the top of the bar
        #         axes[2, i].text(
        #             bar.get_x() + bar.get_width() / 2, 
        #             produced[wished_scenarios].max().max()*.9,
        #             # height + (0.2 * height),  # Adjust the vertical position
        #             label,
        #             ha='center', va='bottom',
        #             fontsize=10, color='black'
        #         )
                
        # y_value =  produced.loc[year, "base"]
        # axes[2,i].axhline(y=y_value, 
        #                   color="black")
        # # Add the text label near the line
        # axes[2, i].text(
        #     x=0.5,  # Position the text at the center of the x-axis (adjust as needed)
        #     y=y_value,  # Align the text vertically with the line
        #     s=f"Base: {y_value:.0f} Mt H$_2$",  # The text to display
        #     color="black",  # Text color
        #     ha='center',  # Horizontal alignment
        #     va='bottom',  # Vertical alignment, place text above the line
        #     fontsize=10  # Adjust the font size as needed
        # )
        
       
       
        
        # # captured CO2 from DAC -----------------------------------------------
        # captured = (dac.unstack()[year].T[["fast", "slow"]]*-1)/1e6
        # captured.plot(kind="bar", legend=False, ax=axes[2, i ],
        #               color=snakemake.config["plotting"]["tech_colors"]["co2"])
        # y_value = ((dac.unstack()[year]).T["base"]/1e6*-1)
        # axes[2,i].axhline(y=y_value, 
        #                   color="black")
        # # Add the text label near the line
        # axes[2, i].text(
        #     x=0.5,  # Position the text at the center of the x-axis (adjust as needed)
        #     y=y_value,  # Align the text vertically with the line
        #     s=f"Base: {y_value:.0f} Mt CO$_2$",  # The text to display
        #     color="black",  # Text color
        #     ha='center',  # Horizontal alignment
        #     va='bottom',  # Vertical alignment, place text above the line
        #     fontsize=10  # Adjust the font size as needed
        # )
        
        # total system costs ------------------------------------------------
        # cost_df.drop(to_drop, errors="ignore", inplace=True)
        diff = (cost_df-cost_df["base"])[wished_scenarios].xs(year, level=1, axis=1)
        diff = diff.rename(index=rename_techs).groupby(level=0).sum()/1e9
        diff = diff[(abs(diff)>1).any(axis=1)]
        
        # Calculate the percentage increase/decrease
        percent_change = cost_df.sum().div(cost_df["base"].sum()).unstack()[year].loc[wished_scenarios]
        
        # diff.T.plot(kind="bar", stacked=True, ax=axes[2,i],
        #           color = [snakemake.config["plotting"]["tech_colors"][i] for i in diff.index],
        #           legend=False
        #           )
        # diff.sum().rename("net").plot(ax=axes[2,i], legend=False)
        # bars = diff.sum().plot(kind="bar", legend=False, ax=axes[3,i], color="gray")
        bars = ((percent_change-1)*100).plot(kind="bar", legend=False, ax=axes[3,i], color="gray")
        
        # Calculate the percentage increase/decrease
        percent_change_labels = [f"{(p-1)*100:.0f}%" for p in percent_change]
        diff_change_labels = [f"{p:.0f}bn Euro" for p in diff.sum()]
        
        # Add text annotations on each bar
        for bar, label in zip(bars.patches, diff_change_labels):
            height = bar.get_height()
            # Position the text slightly above the top of the bar
            axes[3, i].text(
                bar.get_x() + bar.get_width() / 2, 
                height + (0.02 * height),  # Adjust the vertical position
                label,
                ha='center', va='bottom',
                fontsize=10, color='black'
            )
        
    
    # Add row subtitles
    row_titles = [
        "Total Renewable Capacities Difference",
        "Electrified Individual Heating Difference",
        "Hydrogen Production from Electrolysis Difference",
        "Total System Costs Difference"
    ]
    
    for i, title in enumerate(row_titles):
        fig.text(
            0.5,  # Centered on the x-axis
            0.95 - (i * 0.23),  # y-coordinate adjusts for each row, depending on spacing
            title,
            ha='center', va='bottom',
            fontsize=14, # fontweight='bold'
        )
    
    axes[0, 0].set_ylabel("GW")
    axes[0, ncols-1].legend(ncols=1, bbox_to_anchor=(1.8,1), loc="upper right")
    axes[1, 0].set_ylabel("million heat pumps")
    axes[1, ncols-1].legend(ncols=1, bbox_to_anchor=(2.,1), loc="upper right")
    axes[2, 0].set_ylabel("Mt$_{H2}$")
    axes[3, 0].set_ylabel("%")
    
    fig.savefig(snakemake.output.balances[:-19] + "comparison-new.pdf",
                bbox_inches="tight")



def plot_sensi_dsm():
    rename_dict = {"battery": "battery",
                   "home battery": "battery",
                   "rural water tanks": "water tanks",
                   "urban central water tanks": "water tanks",
                   "urban decentral water tanks": "water tanks",
                   }
    stores = ["H2 Store", "battery", "water tanks"]
    generators = ["solar", "solar rooftop", 'solar-hsat', "onwind", "offwind-ac",
                  "offwind-dc"]
    
    caps = capacities.loc["stores"].droplevel([1,2,3],axis=1)[["base", "dsm", "v2g"]]
    to_plot = caps.rename(index=rename_dict).groupby(level=0).sum().loc[stores]
    diff = to_plot.subtract(to_plot.base, axis=0)
    percentage = (to_plot.div(to_plot.base, axis=0)-1)*100
    period = ['2030', '2035', '2040', '2050']
    
    fig, ax = plt.subplots(nrows=1, ncols=len(diff.index),
                           figsize=(12, 2.5), sharex=True)
    for col, store in enumerate(diff.index):
        legend=True if col==2 else False
        (percentage.loc[store].unstack().T.loc[period, ["dsm", "v2g"]]).plot(kind="bar", title=store,
                                                         ax=ax[col],
                                                         legend=legend)
        # write absolute nuumbers on top of the bars
        absolute = diff.loc[store].unstack().T.loc[period, ["dsm", "v2g"]]/1e3
        data_to_plot = percentage.loc[store].unstack().T.loc[period, ["dsm", "v2g"]]
        text = round(absolute).values.flatten(order="K")
        # Add text on top of each bar
        for i, bar in enumerate(ax[col].patches):
            if abs(bar.get_height()) > 2:
                # Add absolute value text
                ax[col].annotate(
                    f'{text[i]} GWh',
                    (bar.get_x() + bar.get_width() / 2, bar.get_height()),  # Position at the top of the bar
                    ha='center', va='bottom', fontsize=9, color='black'  # Styling
                )
            
    ax[0].set_ylabel("difference in energy capacity [%]")
    ax[1].set_xlabel("planning horizon")
    ax[0].set_xlabel("")
    ax[2].set_xlabel("")
    
    fig.savefig(snakemake.output.balances[:-19] + "stores_dsm_v2g.pdf",
                bbox_inches="tight")
    
    
    gens = capacities.loc["generators"].droplevel([1,2,3], axis=1)[["base", "dsm", "v2g"]]
    diff = gens.subtract(gens.base, axis=0)
    to_plot = diff.loc[generators][["dsm", "v2g"]]
    
    fig, ax = plt.subplots(nrows=1, ncols=len(to_plot.index),
                           figsize=(10, 3), sharex=True)
    for col, gen in enumerate(to_plot.index):
        (diff.loc[store].unstack().T.loc[period, ["dsm", "v2g"]]).plot(kind="bar", title=store,
                                                         ax=ax[col])

    ax[0].set_ylabel("difference in energy capacity [%]")
    ax[1].set_xlabel("planning horizon")
    ax[0].set_xlabel("")
    ax[2].set_xlabel("")
    
    fig.savefig(snakemake.output.balances[:-19] + "stores_dsm_v2g.pdf",
                bbox_inches="tight")
    
  
def plot_vehicle_costs():
    
    fn = "/home/lisa/Documents/own_projects/endogenous_transport/data/total_costs.csv"
    new_costs = pd.read_csv(fn, index_col=[0,1])
    
    rename_engine = {"fuel_cell": "FCEV",
                     "ice": "ICE",
                     "electric": "BEV"}
    to_plot = new_costs["fixed"].unstack()
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10,3), sharey=True)
    for i, engine in enumerate(car_keys.keys()):
        a = to_plot[car_keys[engine]]
        (a/1e3).plot(kind="bar", ax=ax[i], title=rename_engine[engine], legend=False,
                     zorder=1)
        # (a/1e3).mean(axis=1).plot(ax=ax[i], kind="bar")
        ax[i].grid(axis="y", zorder=0)
    
    ax[0].set_ylabel("annualised investment costs \n [tEur/a]")
    # Customize the legend
    handles, labels = ax[2].get_legend_handles_labels()
    # Remove "BEV" from each label
    labels = [label.replace("BEV", "").strip() for label in labels]
    ax[2].legend(handles, labels, bbox_to_anchor=(1.05, 0.5))
    fig.savefig(snakemake.output.balances[:-19] + "heavy-duty-annualised-investmentcosts.pdf",
                bbox_inches="tight")
    
    to_plot = new_costs["investment"].unstack()
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10,3), sharey="row")
    for i, engine in enumerate(car_keys.keys()):
        a = to_plot[car_keys[engine]]
        (a/1e3).plot(kind="bar", ax=ax[i], title=rename_engine[engine], legend=False,
                     zorder=1)
        # (a/1e3).mean(axis=1).plot(ax=ax[i], kind="bar")
        ax[i].grid(axis="y", zorder=0)
        # (a.div(a.mean(axis=1), axis=0)-1).plot(kind="bar", ax=ax[1, i], title="deviation mean", legend=False)
        # Add horizontal lines for each investment period's average
        # Scatter the average values for each investment period
        avg_values = (a / 1e3).mean(axis=1)
        for j, period in enumerate(avg_values.index):
            ax[i].scatter(j, avg_values[period], color="black", zorder=3,
                          label="Average" if j == 0 else "")

        
    ax[0].set_ylabel("investment costs \n [tEur]")
    # Customize the legend
    handles, labels = ax[2].get_legend_handles_labels()
    # Remove "BEV" from each label
    labels = [label.replace("BEV", "").strip() for label in labels]
    ax[2].legend(handles, labels, bbox_to_anchor=(1.05, 0.5))
    fig.savefig(snakemake.output.balances[:-19] + "heavy-duty-investmentcosts.pdf",
                bbox_inches="tight")
        
def plot_scenarios(balances, shares):
    carriers = ["land transport demand light", "land transport demand heavy"]
    wished_scenarios = ["low-demand", "base", "high-demand"]
    
    exogen = shares.droplevel([1,2,3], axis=1)
    
    
    fig, axes = plt.subplots(
        nrows=len(carriers), ncols=1, 
        figsize=(14, 7),  sharex=True,
        )
    for row, carrier in enumerate(carriers):
        legend = (row == 0)

        transport_type = carrier.split(" demand")[1]
        exogen_grouped = exogen[exogen.index.str.contains(transport_type)].sum().unstack()
        
        exogen_share = exogen_grouped["existing"].div(exogen_grouped["demand"]).unstack().T
        
        demand = -1 * balances.loc[carrier, "loads"].iloc[0].droplevel([1,2,3]).unstack()
        demand = demand.reindex(wished_scenarios).T/1e9
        
        # plot demand 
        demand.plot(kind="bar", ax=axes[row], title=transport_type[1:],
                    legend=legend, alpha=0.7)
        axes[row].set_ylabel("demand \n [trillion vehicle-km driven]")
        axes[row].set_xlabel("year")
        
        if legend:
            axes[row].legend(loc="upper center", bbox_to_anchor=(0.5, 0.95),
                             ncol=3, frameon=False)
    
        
        # Create secondary y-axis for exogen_share
        ax2 = axes[row].twinx()
        ((1-exogen_share[["fast", "base", "slow"]])*100).plot(ax=ax2,
                                                          style=["--", "-", ":"],
                                                          color="black",
                                                           legend=legend
                                                          )
        # Adjust the legend position
        if legend:  
            ax2.legend(loc="upper left")
        ax2.set_ylabel("Maximum zero emission \n vehicles share \n [%]")
        
        # Optionally, format secondary y-axis label on the right side
        ax2.yaxis.label.set_color("gray")
        
        fig.savefig(snakemake.output.balances[:-19] + "scenario_overview.pdf",
                    bbox_inches="tight")
    
#%%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("plot_summary_all",
                                   configfiles="/home/lisa/Documents/playground/pypsa-eur/config/config.transport_zecm_v2.yaml",)

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    
    n_header = 4
    
    prefix = snakemake.config["run"]["prefix"]
    path = snakemake.output[0].split("graphs")[0]
    scenarios = snakemake.config["run"]["name"]
    
    costs = {}
    balances = {}
    prices = {}
    capacities = {}
    shares = {}
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
            capacities[scenario] = pd.read_csv(f"{path}/{scenario}/csvs/capacities.csv",
                                                      index_col=[0,1],
                                                      header=list(range(n_header)))
            shares[scenario] = pd.read_csv(f"{path}/{scenario}/csvs/shares.csv",
                                                      index_col=[0],
                                                      header=list(range(n_header+1)))
        except FileNotFoundError:
            logger.info(f"{scenario} not solved yet.")
            
    costs = pd.concat(costs, axis=1)
    balances = pd.concat(balances, axis=1)
    prices = pd.concat(prices, axis=1)
    capacities = pd.concat(capacities, axis=1)
    shares = pd.concat(shares, axis=1)
    
    costs.to_csv(snakemake.output.costs_csv)
    balances.to_csv(snakemake.output.balances_csv)
    prices.to_csv(snakemake.output.prices_csv)
      
    plot_costs(costs, drop=True)

    plot_balances(balances, drop=True)
    
    plot_prices(prices, drop=True)
    
    plot_comparison(balances, capacities, drop=True)
    
    # plot_shares_area(balances, shares)


