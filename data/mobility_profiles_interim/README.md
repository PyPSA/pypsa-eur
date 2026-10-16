# Temporary tracked copies of BASt per-segment traffic profiles

**This folder is a temporary stop-gap, not a permanent data source.**

`pkw.csv`, `mot.csv`, `bus.csv`, `lfw.csv`, and `lkw.csv` are weekly traffic-count profiles from the Federal Highway Research Institute (BASt), the same dataset already documented as `mobility_profiles` in `doc/data_inventory.csv` (license: CC-BY-4.0) and tracked in `data/versions.csv`. They are the five per-vehicle-segment profiles consumed by segmented land-transport mode (`sector.transport_segmentation_enable: true`) in `scripts/build_transport_demand.py`.

Of these six, only `mot.csv`, `bus.csv`, `lfw.csv`, and `lkw.csv` are not yet published on the hosted `data.pypsa.org` archive (unlike `kfz.csv`/`pkw.csv`, which already are, via `rules/retrieve.smk`'s `rule retrieve_mobility_profiles`). `pkw.csv` and `kfz.csv` are duplicated here too, purely so this folder could become the archived version to inlcude both lumped mode and segmented mode required files.

For the moment, this folder is where `rules/build_sector.smk`'s `input_transport_traffic_data` reads for segmented mode right now. `kfz.csv` keeps sourcing it directly from the archive via `MOBILITY_PROFILES_DATASET`, exactly as today.

## Removing this folder

Once the hosted archive at `data.pypsa.org/workflows/eur/mobility_profiles/` is extended to publish
`mot.csv`, `bus.csv`, `lfw.csv`, and `lkw.csv` alongside the existing `kfz.csv`/`pkw.csv`:

1. Delete this folder.
2. In `rules/build_sector.smk`, update `input_transport_traffic_data` to source
   `pkw`/`mot`/`bus`/`lfw`/`lkw` from `MOBILITY_PROFILES_DATASET["folder"]` instead of
   `data/mobility_profiles_interim/`.
3. In `rules/retrieve.smk`, update `retrieve_mobility_profiles` to source `pkw`/`mot`/`bus`/`lfw`/`lkw` from `MOBILITY_PROFILES_DATASET["folder"].
4. Remove the `!/data/mobility_profiles_interim/` line from `.gitignore`.
5. Remove the corresponding `[[annotations]]` block from `REUSE.toml`.
