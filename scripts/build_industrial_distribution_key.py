
import pypsa
import pandas as pd
import geopandas as gpd
from shapely import wkt, prepared
from scipy.spatial import cKDTree as KDTree


def prepare_hotmaps_database():

    df = pd.read_csv(snakemake.input.hotmaps_industrial_database,
                     sep=";",
                     index_col=0)

    #remove those sites without valid geometries
    df.drop(df.index[df.geom.isna()],
            inplace=True)

    #parse geometry
    #https://geopandas.org/gallery/create_geopandas_from_pandas.html?highlight=parse#from-wkt-format
    df["Coordinates"] = df.geom.apply(lambda x : wkt.loads(x[x.find(";POINT")+1:]))

    gdf = gpd.GeoDataFrame(df, geometry='Coordinates')

    europe_shape = gpd.read_file(snakemake.input.europe_shape).loc[0, 'geometry']
    europe_shape_prepped = prepared.prep(europe_shape)
    not_in_europe = gdf.index[~gdf.geometry.apply(europe_shape_prepped.contains)]
    print("Removing the following industrial facilities since they are not in European area:")
    print(gdf.loc[not_in_europe])
    gdf.drop(not_in_europe,
             inplace=True)

    country_to_code = {
        'Belgium' : 'BE',
        'Bulgaria' : 'BG',
        'Czech Republic' : 'CZ',
        'Denmark' : 'DK',
        'Germany' : 'DE',
        'Estonia' : 'EE',
        'Ireland' : 'IE',
        'Greece' : 'GR',
        'Spain' : 'ES',
        'France' : 'FR',
        'Croatia' : 'HR',
        'Italy' : 'IT',
        'Cyprus' : 'CY',
        'Latvia' : 'LV',
        'Lithuania' : 'LT',
        'Luxembourg' : 'LU',
        'Hungary' : 'HU',
        'Malta' : 'MA',
        'Netherland' : 'NL',
        'Austria' : 'AT',
        'Poland' : 'PL',
        'Portugal' : 'PT',
        'Romania' : 'RO',
        'Slovenia' : 'SI',
        'Slovakia' : 'SK',
        'Finland' : 'FI',
        'Sweden' : 'SE',
        'United Kingdom' : 'GB',
        'Iceland' : 'IS',
        'Norway' : 'NO',
        'Montenegro' : 'ME',
        'FYR of Macedonia' : 'MK',
        'Albania' : 'AL',
        'Serbia' : 'RS',
        'Turkey' : 'TU',
        'Bosnia and Herzegovina' : 'BA',
        'Switzerland' : 'CH',
        'Liechtenstein' : 'AT',
    }
    gdf["country_code"] = gdf.Country.map(country_to_code)

    if gdf["country_code"].isna().any():
        print("Warning, some countries not assigned an ISO code")

    gdf["x"] = gdf.geometry.x
    gdf["y"] = gdf.geometry.y

    return gdf


def assign_buses(gdf):

    gdf["bus"] = ""

    for c in n.buses.country.unique():
        buses_i = n.buses.index[n.buses.country == c]
        kdtree = KDTree(n.buses.loc[buses_i, ['x','y']].values)

        industry_i = gdf.index[(gdf.country_code == c)]

        if industry_i.empty:
            print("Skipping country with no industry:",c)
        else:
            tree_i = kdtree.query(gdf.loc[industry_i, ['x','y']].values)[1]
            gdf.loc[industry_i, 'bus'] = buses_i[tree_i]

    if (gdf.bus == "").any():
        print("Some industrial facilities have empty buses")
    if gdf.bus.isna().any():
        print("Some industrial facilities have NaN buses")


def build_nodal_distribution_key(gdf):

    sectors = ['Iron and steel','Chemical industry','Cement','Non-metallic mineral products','Glass','Paper and printing','Non-ferrous metals']

    distribution_keys = pd.DataFrame(index=n.buses.index,
                                    columns=sectors,
                                    dtype=float)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout,index_col=0)
    pop_layout["ct"] = pop_layout.index.str[:2]
    ct_total = pop_layout.total.groupby(pop_layout["ct"]).sum()
    pop_layout["ct_total"] = pop_layout["ct"].map(ct_total)
    distribution_keys["population"] = pop_layout["total"]/pop_layout["ct_total"]

    for c in n.buses.country.unique():
        buses = n.buses.index[n.buses.country == c]
        for sector in sectors:
            facilities = gdf.index[(gdf.country_code == c) & (gdf.Subsector == sector)]
            if not facilities.empty:
                emissions = gdf.loc[facilities,"Emissions_ETS_2014"]
                if emissions.sum() == 0:
                    distribution_key = pd.Series(1/len(facilities),
                                                 facilities)
                else:
                    #BEWARE: this is a strong assumption
                    emissions = emissions.fillna(emissions.mean())
                    distribution_key = emissions/emissions.sum()
                distribution_key = distribution_key.groupby(gdf.loc[facilities,"bus"]).sum().reindex(buses,fill_value=0.)
            else:
                distribution_key = distribution_keys.loc[buses,"population"]

            if abs(distribution_key.sum() - 1) > 1e-4:
                print(c,sector,distribution_key)

            distribution_keys.loc[buses,sector] = distribution_key

    distribution_keys.to_csv(snakemake.output.industrial_distribution_key)

if __name__ == "__main__":


    n = pypsa.Network(snakemake.input.network)

    hotmaps_database = prepare_hotmaps_database()

    assign_buses(hotmaps_database)

    build_nodal_distribution_key(hotmaps_database)
