"""Generate scenarios"""


import yaml
import re
from build_energy_totals import build_eea_co2, build_eurostat_co2, build_co2_totals

def emission_sectors_from_opts(opts):

    sectors = ["electricity"]
    if "T" in opts:
        sectors += [
            "rail non-elec",
            "road non-elec"
        ]
    if "H" in opts:
        sectors += [
            "residential non-elec",
            "services non-elec"
        ]
    if "I" in opts:
        sectors += [
            "industrial non-elec",
            "industrial processes",
            "domestic aviation",
            "international aviation",
            "domestic navigation",
            "international navigation"
        ]
    if "A" in opts:
        sectors += [
            "agriculture"
        ]

    return sectors



def co2_emissions_year(countries, opts, year):
    """
    Calculate CO2 emissions in one specific year (e.g. 1990 or 2018).
    """

    eea_co2 = build_eea_co2(year)

    # TODO: read Eurostat data from year > 2014
    # this only affects the estimation of CO2 emissions for BA, RS, AL, ME, MK
    if year > 2014:
        eurostat_co2 = build_eurostat_co2(year=2014)
    else:
        eurostat_co2 = build_eurostat_co2(year)

    co2_totals = build_co2_totals(eea_co2, eurostat_co2)

    sectors = emission_sectors_from_opts(opts)

    co2_emissions = co2_totals.loc[countries, sectors].sum().sum()

    # convert MtCO2 to GtCO2
    co2_emissions *= 0.001

    return co2_emissions





# TODO: move to own rule with sector-opts wildcard?
def build_carbon_budget(o, fn):
    """
    Distribute carbon budget following beta or exponential transition path.
    """
    # opts?

    if "be" in o:
        #beta decay
        carbon_budget = float(o[o.find("cb")+2:o.find("be")])
        be = float(o[o.find("be")+2:])
    if "ex" in o:
        #exponential decay
        carbon_budget = float(o[o.find("cb")+2:o.find("ex")])
        r = float(o[o.find("ex")+2:])

    #TOFIX since no network n
    countries = n.buses.country.dropna().unique()

    e_1990 = co2_emissions_year(countries, opts, year=1990)

    #emissions at the beginning of the path (last year available 2018)
    e_0 = co2_emissions_year(countries, opts, year=2018)

    planning_horizons = snakemake.config['scenario']['planning_horizons']
    t_0 = planning_horizons[0]

    if "be" in o:

        # final year in the path
        t_f = t_0 + (2 * carbon_budget / e_0).round(0)

        def beta_decay(t):
            cdf_term = (t - t_0) / (t_f - t_0)
            return (e_0 / e_1990) * (1 - beta.cdf(cdf_term, be, be))

        #emissions (relative to 1990)
        co2_cap = pd.Series({t: beta_decay(t) for t in planning_horizons}, name=o)

    if "ex" in o:

        T = carbon_budget / e_0
        m = (1 + np.sqrt(1 + r * T)) / T

        def exponential_decay(t):
            return (e_0 / e_1990) * (1 + (m + r) * (t - t_0)) * np.exp(-m * (t - t_0))

        co2_cap = pd.Series({t: exponential_decay(t) for t in planning_horizons}, name=o)

    # TODO log in Snakefile
    if not os.path.exists(fn):
        os.makedirs(fn)
    co2_cap.to_csv(fn, float_format='%.3f')


if __name__ == "__main__":

    print("converting following scenario string to config.yaml:",snakemake.wildcards.sector_opts)

    with open(snakemake.input.config, "r") as input_file:
        config = yaml.safe_load(input_file)

    opts = snakemake.wildcards.sector_opts.split("-")


    config["sector"]["transport"] = "T" in opts
    config["sector"]["heating"] = "H" in opts
    config["sector"]["biomass"] = "B" in opts
    config["sector"]["industry"] = "I" in opts
    config["sector"]["agriculture"] = "A" in opts


    # this case switch should be exhaustive
    for o in opts:
        if o in ["T","H","B","I","A"]:
            pass
        elif o == "nodistrict":
            config["sector"]["district_heating"]["progress"] = 0.0
        elif o == "decentral":
            config["sector"]["electricity_transmission_network"] = False
        elif o[:4] == "wave":
            config["sector"]["wave"] = True
            config["sector"]["wave_cost_factor"] = float(o[4:].replace("p", ".").replace("m", "-"))
            print("Including wave generators with cost factor of", wave_cost_factor)
            add_wave(n, wave_cost_factor)
        elif o[:4] == "dist":
            config["sector"]['electricity_distribution_grid'] = True
            config["sector"]['electricity_distribution_grid_cost_factor'] = float(o[4:].replace("p", ".").replace("m", "-"))
        elif o == "biomasstransport":
                config["sector"]["biomass_transport"] = True
        elif re.match(r'^\d+h$', o, re.IGNORECASE):
            config["snapshots"]["temporal_resolution"] = re.match(r'^\d+h$', o, re.IGNORECASE).group(0)
        elif o[:2] == "cb":
            fn = os.path.join(snakemake.config['results_dir'],
                              snakemake.config['run'],
                              '/csvs/carbon_budget_distribution.csv')
            if not os.path.exists(fn):
                build_carbon_budget(o, fn)
            co2_cap = pd.read_csv(fn, index_col=0).squeeze()
            for year in config["co2_budget"]:
                config["co2_budget"][year] = co2_cap[year]
        elif o[:4] == "Co2L":
            limit = float(o[4:].replace("p", ".").replace("m", "-"))
            config["co2_budget"] = limit
        elif o[:10] == 'linemaxext':
            config["sector"]["line_extension_limit"] = True
            config["sector"]["line_extension_maximum"] = float(o[10:]) * 1e3
        elif "+" in o:
            oo = o.split("+")
            carrier = oo[0]
            attr = oo[1][0]
            #beware if factor is 0 and p_nom_max is np.inf, 0*np.inf is nan
            factor = float(oo[1][1:])

            if carrier not in ["offwind","solar","onwind"]:
                print(f"option {o}: carrier not recognised")
                sys.exit()

            carriers = [carrier]
            if carrier == "offwind":
                carriers = ["offwind-dc","offwind-ac"]

            if attr == "p":
                for c in carriers:
                    config["renewable"][c]["potential_factor"] = factor
            elif attr == "c":
                for c in carriers:
                    config["renewable"][c]["cost_factor"] = factor
            else:
                print(f"option {o}: attribute not recognised")
                sys.exit()
        else:
            print(f"option {o} not recognised")
            sys.exit()

    with open(snakemake.output.config, "w") as output_file:
        yaml.dump(config, output_file)
