import pypsa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pypsa.descriptors import get_switchable_as_dense as as_dense
import geopandas as gpd
from shapely.geometry import Point

solver = "gurobi"

#Load 37 network DE
o = pypsa.Network("./results/networks/elec_s_37_ec_lv1.0_Co2L-3H.nc")
#o = pypsa.examples.scigrid_de(from_master=True)
o.lines.s_max_pu = 0.7
#o.lines.loc[["316", "527", "602"], "s_nom"] = 1715
o.set_snapshots([o.snapshots[12]])

n = o.copy()  # for redispatch model
m = o.copy()  # for market model

o.plot()

#solving original nodal model 
o.lopf(solver_name=solver, pyomo=False)

#######################################
#build market model, with bidding zones 
######################################

#read in shapes
nuts_3_eu = gpd.read_file("data/shapes/NUTS_RG_10M_2021_4326.geojson")
#nuts_1_de = nuts_3_eu.query("LEVL_CODE ==1 and CNTR_CODE == 'DE'")
nuts = nuts_3_eu.query("LEVL_CODE ==1")


#form nodes geometry
#buses_gpd = gpd.GeoDataFrame(n.buses, geometry=gpd.points_from_xy(n.buses.x, n.buses.y),crs="EPSG:4326")
#buses_w_nuts_data = buses_gpd.sjoin(nuts_1_de)

#create one Point for every bus in one column
n.buses['Point'] = [Point(xy) for xy in zip(n.buses.x, n.buses.y)] 

def order_bus_to_state(point):
          for y in nuts['geometry'].tolist():
                  if point.within(y): 
                    value = nuts.loc[ nuts["geometry"] == y, 'NUTS_ID'].iloc[0]
                    return value
                  

#for points in n.buses.items():
#def order_bus_to_state(point):
#          for y in nuts_1_de['geometry'].tolist():
#                  if point.within(y): 
#                    value = nuts_1_de.loc[ nuts_1_de["geometry"] == y, 'NUTS_ID'].iloc[0]
#                    return value
                
                  
states = n.buses["Point"].map(lambda x: order_bus_to_state(x)) 

scen = "BZ5"
def create_scenarios():
    match scen:
            case "BZ2":
                lookup_dict = {"DEF": "DEII1", "DE6": "DEII1", "DE9": "DEII1", "DE3": "DEII1", "DE4": "DEII1",
                          "DE8": "DEII1", "DED": "DEII1", "DEE": "DEII1", "DEG": "DEII1", "DEA": "DEII2",
                          "DEB": "DEII2", "DEC": "DEII2", "DE1": "DEII2", "DE2": "DEII2", "DE7": "DEII2",
                          "OffBZN": "OffBZN", "OffBZB": "OffBZB"}
            case "BZ3":
                lookup_dict = {"DEF": "DEIII1", "DE6": "DEIII1", "DE9": "DEIII1", "DE3": "DEIII2", "DE4": "DEIII2",
                                  "DE8": "DEIII2", "DED": "DEIII2", "DEE": "DEIII2", "DEG": "DEIII2", "DEA": "DEIII3",
                                  "DEB": "DEIII3", "DEC": "DEIII3", "DE1": "DEIII3", "DE2": "DEIII3", "DE7": "DEIII3",
                                  "OffBZN": "OffBZN", "OffBZB": "OffBZB"}
            case "BZ5":
                lookup_dict = {"DEF": "DEV1", "DE6": "DEV2", "DE9": "DEV2", "DE3": "DEV3", "DE4": "DEV3",
                                  "DE8": "DEV3", "DED": "DEV3", "DEE": "DEV3", "DEG": "DEV3", "DEA": "DEV4",
                                  "DEB": "DEV4", "DEC": "DEV4", "DE1": "DEV5", "DE2": "DEV5", "DE7": "DEV5",
                                  "OffBZN": "OffBZN", "OffBZB": "OffBZB"}
    return lookup_dict


zones = states.replace(create_scenarios())
#zones = zones.fillna("?")
#print(zones[zones == "?"])
#wrong = (zones[zones =="?"]).to_frame()
#wrong = wrong.merge(n.buses[['Point']], on="Bus",how='left')
#zones = (n.buses.y > 51).map(lambda x: "North" if x else "South")


for c in m.iterate_components(m.one_port_components):
    c.df.bus = c.df.bus.map(zones)

for c in m.iterate_components(m.branch_components):
    c.df.bus0 = c.df.bus0.map(zones)
    c.df.bus1 = c.df.bus1.map(zones)
    internal = c.df.bus0 == c.df.bus1
    m.mremove(c.name, c.df.loc[internal].index)

m.mremove("Bus", m.buses.index)
m.madd("Bus", zones.unique())
#m.madd("Bus", ["DEII1", "DEII2"])
m.lopf(solver_name=solver, pyomo=False)

m.buses_t.marginal_price

#Build redispatch model n:
p = m.generators_t.p / m.generators.p_nom
n.generators_t.p_min_pu = p
n.generators_t.p_max_pu = p

g_up = n.generators.copy()
g_down = n.generators.copy()

g_up.index = g_up.index.map(lambda x: x + " ramp up")
g_down.index = g_down.index.map(lambda x: x + " ramp down")

up = (
    as_dense(m, "Generator", "p_max_pu") * m.generators.p_nom - m.generators_t.p
).clip(0) / m.generators.p_nom
down = -m.generators_t.p / m.generators.p_nom

up.columns = up.columns.map(lambda x: x + " ramp up")
down.columns = down.columns.map(lambda x: x + " ramp down")

n.madd("Generator", g_up.index, p_max_pu=up, **g_up.drop("p_max_pu", axis=1))

n.madd(
    "Generator",
    g_down.index,
    p_min_pu=down,
    p_max_pu=0,
    **g_down.drop(["p_max_pu", "p_min_pu"], axis=1)
)

#solve redispatch
n.lopf(solver_name=solver, pyomo=False)

#plot both the market results of the x bidding zone market and the redispatch results
fig, axs = plt.subplots(
    1, 3, figsize=(20, 10), subplot_kw={"projection": ccrs.AlbersEqualArea()}
)

market = (
    n.generators_t.p[m.generators.index]
    .T.squeeze()
    .groupby(n.generators.bus)
    .sum()
    .div(2e4)
)
n.plot(ax=axs[0], bus_sizes=market, title="5 bidding zones market simulation")

redispatch_up = (
    n.generators_t.p.filter(like="ramp up")
    .T.squeeze()
    .groupby(n.generators.bus)
    .sum()
    .div(2e4)
)
n.plot(
    ax=axs[1], bus_sizes=redispatch_up, bus_colors="blue", title="Redispatch: ramp up"
)

redispatch_down = (
    n.generators_t.p.filter(like="ramp down")
    .T.squeeze()
    .groupby(n.generators.bus)
    .sum()
    .div(-2e4)
)
n.plot(
    ax=axs[2],
    bus_sizes=redispatch_down,
    bus_colors="red",
    title="Redispatch: ramp down / curtail",
);

#final dispatch looks like:
grouper = n.generators.index.str.split(" ramp", expand=True).get_level_values(0)

n.generators_t.p.groupby(grouper, axis=1).sum().squeeze()


#Higher costs for ramping up/ramping down 
n.generators.loc[n.generators.index.str.contains("ramp up"), "marginal_cost"] *= 2
n.generators.loc[n.generators.index.str.contains("ramp down"), "marginal_cost"] *= -0.5
n.lopf(solver_name=solver, pyomo=False)