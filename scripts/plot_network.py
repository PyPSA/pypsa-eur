
import pandas as pd

#allow plotting without Xwindows
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import pypsa

import numpy as np

from plot_summary import rename_techs, preferred_order

from make_summary import assign_carriers

#from sector/scripts/paper_graphics-co2_sweep.py


from matplotlib.patches import Circle, Ellipse

from matplotlib.legend_handler import HandlerPatch

import cartopy.crs as ccrs


override_component_attrs = pypsa.descriptors.Dict({k : v.copy() for k,v in pypsa.components.component_attrs.items()})
override_component_attrs["Link"].loc["bus2"] = ["string",np.nan,np.nan,"2nd bus","Input (optional)"]
override_component_attrs["Link"].loc["bus3"] = ["string",np.nan,np.nan,"3rd bus","Input (optional)"]
override_component_attrs["Link"].loc["efficiency2"] = ["static or series","per unit",1.,"2nd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["efficiency3"] = ["static or series","per unit",1.,"3rd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["p2"] = ["series","MW",0.,"2nd bus output","Output"]
override_component_attrs["Link"].loc["p3"] = ["series","MW",0.,"3rd bus output","Output"]
override_component_attrs["StorageUnit"].loc["p_dispatch"] = ["series","MW",0.,"Storage discharging.","Output"]
override_component_attrs["StorageUnit"].loc["p_store"] = ["series","MW",0.,"Storage charging.","Output"]



def rename_techs_tyndp(tech):
    tech = rename_techs(tech)
    if "heat pump" in tech or "resistive heater" in tech:
        return "power-to-heat"
    elif tech in ["methanation","hydrogen storage","helmeth"]:
        return "power-to-gas"
    elif tech in ["OCGT","CHP","gas boiler"]:
        return "gas-to-power/heat"
    elif "solar" in tech:
        return "solar"
    elif tech == "Fischer-Tropsch":
        return "power-to-liquid"
    elif "offshore wind" in tech:
        return "offshore wind"
    else:
        return tech

def make_handler_map_to_scale_circles_as_in(ax, dont_resize_actively=False):
    fig = ax.get_figure()
    def axes2pt():
        return np.diff(ax.transData.transform([(0,0), (1,1)]), axis=0)[0] * (72./fig.dpi)

    ellipses = []
    if not dont_resize_actively:
        def update_width_height(event):
            dist = axes2pt()
            for e, radius in ellipses: e.width, e.height = 2. * radius * dist
        fig.canvas.mpl_connect('resize_event', update_width_height)
        ax.callbacks.connect('xlim_changed', update_width_height)
        ax.callbacks.connect('ylim_changed', update_width_height)

    def legend_circle_handler(legend, orig_handle, xdescent, ydescent,
                              width, height, fontsize):
        w, h = 2. * orig_handle.get_radius() * axes2pt()
        e = Ellipse(xy=(0.5*width-0.5*xdescent, 0.5*height-0.5*ydescent), width=w, height=w)
        ellipses.append((e, orig_handle.get_radius()))
        return e
    return {Circle: HandlerPatch(patch_func=legend_circle_handler)}



def make_legend_circles_for(sizes, scale=1.0, **kw):
    return [Circle((0,0), radius=(s/scale)**0.5, **kw) for s in sizes]

def assign_location(n):
    for c in n.iterate_components(n.one_port_components|n.branch_components):

        ifind = pd.Series(c.df.index.str.find(" ",start=4),c.df.index)

        for i in ifind.value_counts().index:
            #these have already been assigned defaults
            if i == -1:
                continue

            names = ifind.index[ifind == i]

            c.df.loc[names,'location'] = names.str[:i]



def plot_map(components=["links","stores","storage_units","generators"],bus_size_factor=1.7e10):

    n = pypsa.Network(snakemake.input.network,
                      override_component_attrs=override_component_attrs)


    assign_location(n)

    #Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"],inplace=True)

    costs = pd.DataFrame(index=n.buses.index)

    for comp in components:
        getattr(n,comp)["nice_group"] = getattr(n,comp).carrier.map(rename_techs_tyndp)

        attr = "e_nom_opt" if comp == "stores" else "p_nom_opt"

        costs = pd.concat((costs,(getattr(n,comp).capital_cost*getattr(n,comp)[attr]).groupby((getattr(n,comp).location,getattr(n,comp).nice_group)).sum().unstack().fillna(0.)),axis=1)

        print(comp,costs)
    costs = costs.groupby(costs.columns,axis=1).sum()


    costs.drop(list(costs.columns[(costs == 0.).all()]),axis=1,inplace=True)

    new_columns = (preferred_order&costs.columns).append(costs.columns.difference(preferred_order))

    costs = costs[new_columns]

    #print(costs)
    #print(costs.sum())

    costs = costs.stack()#.sort_index()

    #print(costs)

    fig, ax = plt.subplots(subplot_kw={"projection":ccrs.PlateCarree()})

    fig.set_size_inches(7,6)

    linewidth_factor=2e3
    ac_color="gray"
    dc_color="m"

    #hack because impossible to drop buses...
    n.buses.loc["EU gas",["x","y"]] = n.buses.loc["DE0 0",["x","y"]]

    n.links.drop(n.links.index[(n.links.carrier != "DC") & (n.links.carrier != "B2B")],inplace=True)


    if snakemake.wildcards.lv == "1.0":
        #should be zero
        line_widths_exp = pd.concat(dict(Line=(n.lines.s_nom_opt-n.lines.s_nom), Link=(n.links.p_nom_opt-n.links.p_nom)))
    else:
        line_widths_exp = pd.concat(dict(Line=(n.lines.s_nom_opt-n.lines.s_nom_min), Link=(n.links.p_nom_opt-n.links.p_nom_min)))

    #PDF has minimum width, so set these to zero
    line_threshold = 500.

    line_widths_exp[line_widths_exp < line_threshold] = 0.

    #drop non-bus
    to_drop = costs.index.levels[0]^n.buses.index
    print("dropping non-buses",to_drop)
    costs.drop(to_drop,level=0,inplace=True,axis=0)

    #make sure they are removed from index
    costs.index = pd.MultiIndex.from_tuples(costs.index.values)

    n.plot(bus_sizes=costs/bus_size_factor,
           bus_colors=snakemake.config['plotting']['tech_colors'],
           line_colors=dict(Line=ac_color, Link=dc_color),
           line_widths=line_widths_exp/linewidth_factor,
           ax=ax)



    handles = make_legend_circles_for([5e9, 1e9], scale=bus_size_factor, facecolor="gray")
    labels = ["{} bEUR/a".format(s) for s in (5, 1)]
    l2 = ax.legend(handles, labels,
                       loc="upper left", bbox_to_anchor=(0.01, 1.01),
                       labelspacing=1.0,
                       framealpha=1.,
                       title='System cost',
                       handler_map=make_handler_map_to_scale_circles_as_in(ax))
    ax.add_artist(l2)

    handles = []
    labels = []

    for s in (10, 5):
        handles.append(plt.Line2D([0],[0],color=ac_color,
                                  linewidth=s*1e3/linewidth_factor))
        labels.append("{} GW".format(s))
    l1 = l1_1 = ax.legend(handles, labels,
                              loc="upper left", bbox_to_anchor=(0.30, 1.01),
                              framealpha=1,
                              labelspacing=0.8, handletextpad=1.5,
                              title='Transmission reinforcement')
    ax.add_artist(l1_1)


    #ax.set_title("Scenario {} with {} transmission".format(snakemake.config['plotting']['scenario_names'][flex],"optimal" if line_limit == "opt" else "no"))


    fig.tight_layout()

    fig.savefig(snakemake.output.map,transparent=True)



def plot_h2_map():

    n = pypsa.Network(snakemake.input.network,
                      override_component_attrs=override_component_attrs)

    if "H2 pipeline" not in n.links.carrier.unique():
        return

    assign_location(n)

    bus_size_factor=1e5
    linewidth_factor=1e4
    #MW below which not drawn
    line_threshold=1e3
    bus_color="m"
    link_color="c"

    #Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"],inplace=True)

    elec = n.links.index[n.links.carrier == "H2 Electrolysis"]

    bus_sizes = pd.Series(0., index=n.buses.index)
    bus_sizes.loc[elec.str.replace(" H2 Electrolysis","")] = n.links.loc[elec,"p_nom_opt"].values/bus_size_factor

    #make a fake MultiIndex so that area is correct for legend
    bus_sizes.index = pd.MultiIndex.from_product([bus_sizes.index,["electrolysis"]])

    n.links.drop(n.links.index[n.links.carrier != "H2 pipeline"],inplace=True)

    link_widths = n.links.p_nom_opt/linewidth_factor
    link_widths[n.links.p_nom_opt < line_threshold] = 0.

    n.links.bus0 = n.links.bus0.str.replace(" H2","")
    n.links.bus1 = n.links.bus1.str.replace(" H2","")

    print(link_widths.sort_values())

    print(n.links[["bus0","bus1"]])

    fig, ax = plt.subplots(subplot_kw={"projection":ccrs.PlateCarree()})

    fig.set_size_inches(7,6)

    n.plot(bus_sizes=bus_sizes,
           bus_colors={"electrolysis" : bus_color},
           line_colors=dict(Link=link_color),
           line_widths={"Link" : link_widths},
           branch_components=["Link"],
           ax=ax)

    handles = make_legend_circles_for([50000, 10000], scale=bus_size_factor, facecolor=bus_color)
    labels = ["{} GW".format(s) for s in (50, 10)]
    l2 = ax.legend(handles, labels,
                       loc="upper left", bbox_to_anchor=(0.01, 1.01),
                       labelspacing=1.0,
                       framealpha=1.,
                       title='Electrolyzer capacity',
                       handler_map=make_handler_map_to_scale_circles_as_in(ax))
    ax.add_artist(l2)

    handles = []
    labels = []

    for s in (50, 10):
        handles.append(plt.Line2D([0],[0],color=link_color,
                                  linewidth=s*1e3/linewidth_factor))
        labels.append("{} GW".format(s))
    l1 = l1_1 = ax.legend(handles, labels,
                              loc="upper left", bbox_to_anchor=(0.30, 1.01),
                              framealpha=1,
                              labelspacing=0.8, handletextpad=1.5,
                              title='H2 pipeline capacity')
    ax.add_artist(l1_1)


    #ax.set_title("Scenario {} with {} transmission".format(snakemake.config['plotting']['scenario_names'][flex],"optimal" if line_limit == "opt" else "no"))

    fig.tight_layout()

    fig.savefig(snakemake.output.map.replace("-costs-all","-h2_network"),transparent=True)


def plot_map_without():

    n = pypsa.Network(snakemake.input.network,
                      override_component_attrs=override_component_attrs)

    assign_location(n)

    #Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"],inplace=True)

    fig, ax = plt.subplots(subplot_kw={"projection":ccrs.PlateCarree()})

    fig.set_size_inches(7,6)

    linewidth_factor=2e3
    ac_color="gray"
    dc_color="m"

    #hack because impossible to drop buses...
    n.buses.loc["EU gas",["x","y"]] = n.buses.loc["DE0 0",["x","y"]]

    n.links.drop(n.links.index[(n.links.carrier != "DC") & (n.links.carrier != "B2B")],inplace=True)


    if snakemake.wildcards.lv == "1.0":
        line_widths_exp = pd.concat(dict(Line=(n.lines.s_nom), Link=(n.links.p_nom)))
    else:
        line_widths_exp = pd.concat(dict(Line=(n.lines.s_nom_min), Link=(n.links.p_nom_min)))

    #PDF has minimum width, so set these to zero
    line_threshold = 0.

    line_widths_exp[line_widths_exp < line_threshold] = 0.

    line_widths_exp[line_widths_exp > 1e4] = 1e4


    n.plot(bus_sizes=10,
           bus_colors="k",
           line_colors=dict(Line=ac_color, Link=dc_color),
           line_widths=line_widths_exp/linewidth_factor,
           ax=ax)

    handles = []
    labels = []

    for s in (10, 5):
        handles.append(plt.Line2D([0],[0],color=ac_color,
                                  linewidth=s*1e3/linewidth_factor))
        labels.append("{} GW".format(s))
    l1 = l1_1 = ax.legend(handles, labels,
                              loc="upper left", bbox_to_anchor=(0.05, 1.01),
                              framealpha=1,
                              labelspacing=0.8, handletextpad=1.5,
                              title='Today\'s transmission')
    ax.add_artist(l1_1)


    #ax.set_title("Scenario {} with {} transmission".format(snakemake.config['plotting']['scenario_names'][flex],"optimal" if line_limit == "opt" else "no"))


    fig.tight_layout()

    fig.savefig(snakemake.output.today,transparent=True)


def plot_series(carrier="AC"):

    n = pypsa.Network(snakemake.input.network,
                      override_component_attrs=override_component_attrs)


    assign_location(n)

    assign_carriers(n)

    buses = n.buses.index[n.buses.carrier == carrier]

    supply = pd.DataFrame(index=n.snapshots)
    for c in n.iterate_components(n.branch_components):
        for i in range(2):
            supply = pd.concat((supply,(-1)*c.pnl["p"+str(i)].loc[:,c.df.index[c.df["bus" + str(i)].isin(buses)]].groupby(c.df.carrier,axis=1).sum()),axis=1)

    for c in n.iterate_components(n.one_port_components):
        comps = c.df.index[c.df.bus.isin(buses)]
        supply = pd.concat((supply,((c.pnl["p"].loc[:,comps]).multiply(c.df.loc[comps,"sign"])).groupby(c.df.carrier,axis=1).sum()),axis=1)

    supply = supply.groupby(rename_techs_tyndp,axis=1).sum()

    both = supply.columns[(supply < 0.).any() & (supply > 0.).any()]

    positive_supply = supply[both]
    negative_supply = supply[both]

    positive_supply[positive_supply < 0.] = 0.
    negative_supply[negative_supply > 0.] = 0.

    supply[both] = positive_supply

    suffix = " charging"

    negative_supply.columns = negative_supply.columns + suffix

    supply = pd.concat((supply,negative_supply),axis=1)


    fig, ax = plt.subplots()

    fig.set_size_inches((8,5))

    #14-21.2 for flaute
    #19-26.1 for flaute

    start = "2013-01-19"

    stop = "2013-01-26"

    threshold = 10e3

    to_drop = supply.columns[(abs(supply) < threshold).all()]

    print("dropping",to_drop)

    supply.drop(columns=to_drop,inplace=True)

    supply.index.name=None

    supply = supply/1e3

    supply.rename(columns={"electricity" : "electric demand",
                           "heat" : "heat demand"},
                  inplace=True)

    preferred_order = pd.Index(["electric demand","transmission lines","hydroelectricity","hydro reservoir","run of river","pumped hydro storage","CHP","onshore wind","offshore wind","solar PV","solar thermal","building retrofitting","ground heat pump","air heat pump","resistive heater","OCGT","gas boiler","gas","natural gas","methanation","hydrogen storage","battery storage","hot water storage"])

    new_columns = (preferred_order&supply.columns).append(supply.columns.difference(preferred_order))


    supply.loc[start:stop,new_columns].plot(ax=ax,kind="area",stacked=True,linewidth=0.,color=[snakemake.config['plotting']['tech_colors'][i.replace(suffix,"")] for i in new_columns])



    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    new_handles = []
    new_labels = []

    for i,item in enumerate(labels):
        if "charging" not in item:
            new_handles.append(handles[i])
            new_labels.append(labels[i])

    ax.legend(new_handles,new_labels,ncol=3,loc="upper left")

    ax.set_xlim([start,stop])

    ax.set_ylim([-1300,1900])

    ax.grid(True)

    ax.set_ylabel("Power [GW]")

    fig.tight_layout()

    fig.savefig("{}/series-{}-{}-{}.pdf".format(snakemake.config['summary_dir'],carrier,start,stop),transparent=True)


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
        snakemake.input.scenario = "lv1.0"  #lv1.0, lv1.25, lvopt
        snakemake.config["run"] = "190503-es2050-lv"
        snakemake.input.network = "{}{}/postnetworks/elec_s_181_{}__Co2L0-3H-T-H-B-I-solar3.nc".format(snakemake.config['results_dir'],
                                                                                                       snakemake.config['run'],
                                                                                                       snakemake.input.scenario)
        snakemake.output.network = "{}{}/maps/elec_s_181_{}__Co2L0-3H-T-H-B-I-solar3.pdf".format(snakemake.config['results_dir'],
                                                                                                 snakemake.config['run'],
                                                                                                 snakemake.input.scenario)


    plot_map(components=["generators","links","stores","storage_units"],bus_size_factor=1.5e10)

    plot_h2_map()

    plot_map_without()

    #plot_map(components=["generators"],bus_size_factor=1.7e10,suffix="generators")

    #plot_map(components=["links","stores","storage_units"],bus_size_factor=4e9,suffix="rest")

    #plot_series(carrier="AC")

    #plot_series(carrier="heat")
