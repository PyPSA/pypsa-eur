# supply_curve_v5 — Run Notes

- Used the new DEA cost base for CDR
- Costs are too low and way too much CDR gets deployed, especially DACs

## Config changes vs v4

- `solid_biomass_potential_factor`: raised from `0.75` → `1.0` (full EU biomass potential, no longer capped)
- Additional BECCS pathways enabled: `biogas_upgrading_cc`, `biomass_to_liquid_cc`, `biosng_cc`, `bioH2`, `methanol.biomass_to_methanol_cc`
- Removed erroneous `costs.overwrites` for `biomass CHP capture`: v4 had set `compression-electricity-input = 0.16`, which is actually the DEA value for `compression-heat-output` (wrong field). The correct DEA `compression-electricity-input` is 0.085 MWh/tCO₂. Dropping the overwrite restores the correct value. The `electricity-input = 0.025` overwrite was also removed as it matched the DEA default anyway.
