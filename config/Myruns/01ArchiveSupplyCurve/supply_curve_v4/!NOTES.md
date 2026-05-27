# supply_curve_v4 — Run Notes

- Used the new DEA cost base for CDR
- Costs are too low and way too much CDR gets deployed, especially DACs

## Config changes vs v3

- CDR credit price range €0–€500 in €50 steps (same as v3)
- Temporal resolution: 168 segments (same as v3)
- DAC weather factors enabled (same as v3)
- Biomass transport enabled (new in v4: `sector.biomass_transport: true`)
- CO2 prices updated: 2030=119 (was 140), 2040=279 (was 180), 2050=463 (was 460) €/tCO₂
- Gurobi: added `DualReductions: 0` and `InfUnbdInfo: 1` for better solver diagnostics
