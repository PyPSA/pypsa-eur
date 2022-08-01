##########################################
Limitations
##########################################

While the benefit of an openly available, functional and partially validated
model of the European energy system is high, many approximations have
been made due to missing data.
The limitations of the dataset are listed below,
both as a warning to the user and as an encouragement to assist in
improving the approximations.

This list of limitations is incomplete and will be added to over time.

See also the `GitHub repository issues <https://github.com/PyPSA/pypsa-eur-sec/issues>`_.

- **Electricity transmission network topology:**
  The grid data is based on a map of the ENTSO-E area that is known
  to contain small distortions to improve readability. Since the exact impedances
  of the lines are unknown, approximations based on line lengths and standard
  line parameters were made that ignore specific conductoring choices for
  particular lines. There is no openly available data on busbar configurations, switch
  locations, transformers or reactive power compensation assets.

- **Assignment of electricity demand to transmission nodes:**
  Using Voronoi cells to aggregate load and generator data to transmission
  network substations ignores the topology of the underlying distribution network,
  meaning that assets may be connected to the wrong substation.

- **Incomplete information on existing assets:** Approximations have
  been made for missing data, including: existing distribution grid
  capacities and costs, existing space and water heating supply,
  existing industry facilities, existing transport vehicle fleets.

- **Exogenous pathways for transformation of transport and industry:**
  To avoid penny-switching the transformation of transport and
  industry away from fossil fuels is determined exogenously.

- **Energy demand distribution within countries:**
  Assumptions
  have been made about the distribution of demand in each country proportional to
  population and GDP that may not reflect local circumstances.
  Openly available
  data on load time series may not correspond to the true vertical load and is
  not spatially disaggregated; assuming, as we have done, that the load time series
  shape is the same at each node within each country ignores local differences.

- **Hydro-electric power plants:**
  The database of hydro-electric power plants does not include plant-specific
  energy storage information, so that blanket values based on country storage
  totals have been used. Inflow time series are based on country-wide approximations,
  ignoring local topography and basin drainage; in principle a full
  hydrological model should be used.

- **International interactions:**
  Border connections and power flows to Russia,
  Belarus, Ukraine, Turkey and Morocco have not been taken into account;
  islands which are not connected to the main European system, such as Malta,
  Crete and Cyprus, are also excluded from the model.

- **Demand sufficiency:** Further measures of demand reduction may be
  possible beyond the assumptions made here.
