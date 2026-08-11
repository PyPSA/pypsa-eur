<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Weather Data {#cutouts}

Cutouts are spatio-temporal subsets of the European weather data from the [ECMWF ERA5](https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation) reanalysis dataset and the [CMSAF SARAH-3](https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002) solar surface radiation dataset.
They have been prepared by and are for use with the [atlite](https://github.com/PyPSA/atlite) tool.
The [tutorial](tutorial.md) uses a smaller cutout than required for the full model (30 MB), which is also automatically downloaded.

There are two ways to obtain cutouts:

- **Retrieve pre-built cutouts** from the archive using the `retrieve_cutout` rule (default). This downloads ready-made cutouts from `data.pypsa.org`.
- **Build cutouts from scratch** using the `build_cutout` rule. This requires access to the [CDS API](https://cds.climate.copernicus.eu/api-how-to) to download ERA5 data directly.

!!! info "See also"
    For building your own cutouts, see [build_cutout][] and the [atlite documentation](https://atlite.readthedocs.io).

The following pre-built cutouts are available for download under
`https://data.pypsa.org/workflows/cutout/<version>/<cutout>.nc`
(click on a cutout name below to download directly).

**v1.0**

| Cutout | Size |
|--------|------|
| [be-03-2013-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/be-03-2013-era5.nc) | 11.1 MB |
| [dach-03-2013-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/dach-03-2013-sarah3-era5.nc) | 35.1 MB |
| [europe-1995-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-1995-sarah3-era5.nc) | 6.2 GB |
| [europe-1996-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-1996-sarah3-era5.nc) | 6.5 GB |
| [europe-2008-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2008-sarah3-era5.nc) | 6.2 GB |
| [europe-2009-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2009-sarah3-era5.nc) | 6.2 GB |
| [europe-2010-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2010-sarah3-era5.nc) | 6.1 GB |
| [europe-2012-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2012-sarah3-era5.nc) | 6.6 GB |
| [europe-2013-03-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2013-03-sarah3-era5.nc) | 140.5 MB |
| [europe-2013-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2013-sarah3-era5.nc) | 6.1 GB |
| [europe-2019-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2019-sarah3-era5.nc) | 6.1 GB |
| [europe-2020-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2020-sarah3-era5.nc) | 6.6 GB |
| [europe-2021-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2021-sarah3-era5.nc) | 6.1 GB |
| [europe-2023-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2023-sarah3-era5.nc) | 6.1 GB |
| [europe-2024-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2024-sarah3-era5.nc) | 6.7 GB |
| [europe-2025-sarah3-era5.nc](https://data.pypsa.org/workflows/cutout/v1.0/europe-2025-sarah3-era5.nc) | 6.7 GB |

**Relevant Settings**

```yaml
atlite:
  default_cutout:
  cutouts:
```

!!! info "See also"
    Documentation of the configuration file `config/config.yaml` at
    [atlite](configuration.md#atlite_cf) and [data](configuration.md#data_cf).

**Outputs**

- `cutouts/{cutout}`: weather data from either the [ERA5](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) reanalysis weather dataset and/or [SARAH-3](https://wui.cmsaf.eu/safira/action/viewProduktSearch) satellite-based historic weather data.
