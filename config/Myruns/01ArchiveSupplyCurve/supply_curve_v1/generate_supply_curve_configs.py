"""
Generate supply curve scenario configs for CDR credit prices €0–€750 in €75 steps.

Each priced config (S01-S10) is identical to S0 except:
- cdr_credit_scope, cdr_credit_limit_by_year, and cdr_credit_prices_by_scope are added
- Credit price is flat across all years (same value for 2025/2030/2040/2050)
- Credit limit follows recommended market demand caps:
    2030: 100 Mt/yr (consistent with NZIA storage ramp-up)
    2040: 300 Mt/yr (EC ICM Strategy anchor)
    2050: 500 Mt/yr (CATF upper-bound storage capacity)
- Planning horizons: 2030, 2040, 2050

S0 remains a special case with `cdr_credit_timing: capture` so the zero-price
anchor does not activate solve-time provenance tracking.

Usage:
    python generate_supply_curve_configs.py
"""

import os

TEMPLATE = """\
# ==========================================
# S{index:02d} — CDR CREDIT PRICE €{price}/tCO₂ (flat across all years)
# Supply curve scenario: credit price = €{price}/tCO₂ for DAC and BECCS.
# Credit demand cap: 100 Mt/yr (2030), 300 Mt/yr (2040), 500 Mt/yr (2050)
# based on NZIA / EC ICM Strategy / CATF recommended market demand levels.
# Planning horizons: 2030, 2040, 2050.
# Each horizon uses its own previous solve as brownfield (2030→2040→2050),
# so CDR capacity built in 2030 is correctly inherited into 2040.
# ==========================================

run:
  name: "S{index:02d}-cdr-{price:03d}eur-336seg"
  disable_progressbar: true
  use_shadow_directory: false
  shared_resources:
    policy: false

foresight: myopic

scenario:
  clusters: [96]
  opts: ['']
  sector_opts: ['336seg']
  planning_horizons: [2030, 2040, 2050]

countries:
  ['AL','AT','BA','BE','BG','CH','CZ','DE','DK','EE','ES','FI','FR','GB','GR',
   'HR','HU','IE','IT','LT','LU','LV','ME','MD','MK','NL','NO','PL','PT',
   'RO','RS','SE','SI','SK','UA','XK']

snapshots:
  start: "2019-01-01"
  end: "2020-01-01"
  inclusive: "left"

atlite:
  default_cutout: "europe-2019-era5"
  cutouts:
    "europe-2019-era5":
      module: era5
      x: [-12.0, 42.0]
      y: [33.0, 72.0]
      dx: 0.3
      dy: 0.3
      time: ["2019", "2019"]
      chunks:
        time: 500
      prepare_kwargs:
        monthly_requests: true
        tmpdir: ./cutouts_tmp/

electricity:
  co2limit_enable: false
  transmission_limit: c1.5

  extendable_carriers:
    Generator: [nuclear, solar, solar-hsat, onwind, offwind-ac, offwind-dc, offwind-float, OCGT, CCGT, biomass]
    StorageUnit: [battery]
    Store: [H2]
    Link: []

  conventional_carriers: [oil, OCGT, CCGT, coal, lignite, geothermal, biomass]
  renewable_carriers: [solar, solar-hsat, onwind, offwind-ac, offwind-dc, offwind-float, hydro]

  estimate_renewable_capacities:
    enable: true
    from_gem: true
    year: 2020
    expansion_limit: false
    technology_mapping:
      Offshore: offwind-ac
      Onshore: onwind
      PV: solar

conventional:
  unit_commitment: false
  dynamic_fuel_price: true
  fuel_price_rolling_window: 6
  nuclear:
    p_max_pu: data/nuclear_p_max_pu.csv
  biomass:
    p_max_pu: 0.65

sector:
  biomass: true
  industry: true
  dac: true
  dac_variants:
    enable: true
    variants:
      liquid:
        carrier_name: DAC-liquidHT
        cost_technology: direct air capture - liquid
      solid:
        carrier_name: DAC-solidLT
        cost_technology: direct air capture - solid
      electrochemical:
        carrier_name: DAC-electrochemical
        cost_technology: direct air capture - electrochemical
  cdr_credit_scope: ["dac", "biogenic"]

  cdr_credit_limit_by_year:   # market demand caps (Mt CO₂/yr)
    2030: 100   # NZIA storage ramp-up
    2040: 300   # EC ICM Strategy anchor
    2050: 500   # CATF upper-bound storage capacity

  cdr_credit_prices_by_scope:
    dac:
      2025: {price}
      2030: {price}
      2040: {price}
      2050: {price}
    biogenic:
      2025: {price}
      2030: {price}
      2040: {price}
      2050: {price}

  co2_vent: false
  co2_spatial: true
  co2_network: true
  co2_network_cost_factor:
    2030: 5.0
    2040: 3.0
    2050: 1.5
  cc_fraction: 0.9

  # Credits are issued at geological sequestration, not at capture
  # This prevents overcrediting CO2 that is later re-emitted via e-fuels
  cdr_credit_timing: sequestration

  # Standalone CDR market: DAC/BECCS earn only CDR credit, not CO2 price on top
  cdr_credit_standalone: true

  # EU-ambition CO2 T&S infrastructure ramp-up (Mt CO₂/yr)
  # Consistent with EC CCUS strategy: ~50 Mt by 2030, scaling to ~550 Mt by 2050
  co2_sequestration_potential:
    2025: 50
    2030: 150
    2035: 275
    2040: 400
    2045: 700
    2050: 1000

  # Use full ENSPRESO solid biomass potential (~1020 TWh/yr by 2040)
  # Restored after reduced 50% cap made 2040 brownfield solves infeasible
  solid_biomass_potential_factor: 0.75
  # CCS deployment rate limits — prevent unrealistic overnight build-up.
  # Calibrated against EU CCS project pipelines and published transition pathways:
  #   process emissions CC : ~20 MtCO2/yr by 2030, growing to ~60 by 2050
  #   SMR CC               : ~8 MtCO2/yr by 2030 (blue H2); falls naturally
  #                          with higher CDR credit prices
  #   urban central gas    : capped to prevent new fossil gas CHP lock-in in 2050;
  #   CHP CC                 reduces from current ~96 MtCO2/yr to ≤15 MtCO2/yr
  # factor 1.3 is the built-in 30% slack buffer in add_max_growth().
  # Renewable limits are carried over from config.default.yaml defaults.
  limit_max_growth:
    enable: true
    factor: 1.3
    max_growth:
      onwind: 16
      solar: 28
      "solar-hsat": 7
      "solar rooftop": 15
      "offwind-ac": 35
      "offwind-dc": 35
      "process emissions CC": 0.18
      "SMR CC": 0.39
      "urban central gas CHP CC": 0.25
    max_relative_growth:
      onwind: 3
      solar: 3
      "solar-hsat": 3
      "solar rooftop": 3
      "offwind-ac": 3
      "offwind-dc": 3

  allam_cycle_gas: false
  enhanced_geothermal:
    enable: false

  regional_co2_sequestration_potential:

    enable: true
    attribute:
    - conservative estimate Mt
    - conservative estimate GAS Mt
    - conservative estimate OIL Mt
    - conservative estimate aquifer Mt
    include_onshore: false
    min_size: 3
    max_size: 25
    years_of_storage: 25

  co2_sequestration_cost: 20
  co2_sequestration_cost_by_country:
    "AL": 12
    "BE": 7
    "BG": 13
    "DE": 8
    "DK": 6
    "EE": 13
    "ES": 10
    "FI": 12
    "FR": 9
    "GB": 6
    "HR": 11
    "IE": 8
    "IT": 10
    "LT": 13
    "LV": 13
    "ME": 12
    "NL": 6
    "NO": 5
    "PL": 11
    "PT": 11
    "RO": 12
    "SE": 10
    "UA": 12
  co2_sequestration_lifetime: 50

co2_budget: null  # No hard CO2 cap — CDR deploys on market economics, not mandate

clustering:
  temporal:
    resolution_elec: false
    resolution_sector: "336seg"

lines:
  reconnect_crimea: true
  dynamic_line_rating:
    activate: false

costs:
  overwrites:
    compression-electricity-input:
      biomass CHP capture: 0.16
    electricity-input:
      biomass CHP capture: 0.025
  emission_prices:
    enable: true
    dynamic: false
    co2:
      2030: 140
      2040: 180
      2050: 460

solving:
  options:
    io_api: mps
  agg_p_nom_limits:
    agg_offwind: false
    agg_solar: false
    include_existing: false
    file: config/Myruns/supply_curve/nuclear_p_nom_limits.csv
  constraints:
    CCL: true
    EQ: false
    BAU: false
    SAFE: false
  solver:
    name: gurobi
    options: gurobi-default

  solver_options:
    gurobi-default:
      threads: 8
      method: 2
      crossover: 0
      NumericFocus: 3
      BarHomogeneous: 1
      ScaleFlag: 2
      ObjScale: -0.5
      BarConvTol: 1.0e-5
      FeasibilityTol: 1.0e-5
      OptimalityTol: 1.0e-5
      Seed: 123
      GURO_PAR_BARDENSETHRESH: 200
      TimeLimit: 172800

data:
  cutout:
    source: build
    version: unknown
  wdpa:
    source: archive
    version: latest
  wdpa_marine:
    source: archive
    version: latest
"""

output_dir = os.path.dirname(os.path.abspath(__file__))
prices = range(75, 751, 75)  # 75, 150, 225, ..., 750 (S0 at €0 already exists)


def ensure_s0_capture_timing(filepath: str) -> bool:
    """Keep S0 aligned with the zero-price shortcut used in the thesis runs."""
    with open(filepath, encoding="utf-8") as f:
        content = f.read()

    capture_line = "  cdr_credit_timing: capture"
    sequestration_line = "  cdr_credit_timing: sequestration"

    if capture_line in content:
        return False
    if sequestration_line not in content:
        raise ValueError(f"Could not find cdr_credit_timing setting in {filepath}")

    updated = content.replace(sequestration_line, capture_line, 1)
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(updated)
    return True


s0_path = os.path.join(output_dir, "config.S0.yaml")
s0_updated = ensure_s0_capture_timing(s0_path)

for index, price in enumerate(prices, start=1):
    filename = f"config.S{index:02d}-{price:03d}eur.yaml"
    filepath = os.path.join(output_dir, filename)
    content = TEMPLATE.format(index=index, price=price)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Created {filename}")

print(f"\nDone. {len(list(prices))} scenario configs created (S01–S{len(list(prices)):02d}).")
print("S00 (€0) = config.S0.yaml maintained with capture-time crediting.")
if s0_updated:
    print("Updated config.S0.yaml to keep cdr_credit_timing: capture.")
print("\nScenario summary:")
print(f"  S00: €0/t   (config.S0.yaml)")
for index, price in enumerate(prices, start=1):
    print(f"  S{index:02d}: €{price}/t  (config.S{index:02d}-{price:03d}eur.yaml)")
