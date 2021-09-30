.. _spatial_resolution:

##########################################
Spatial resolution
##########################################

The default nodal resolution of the model follows the electricity
generation and transmission model `PyPSA-Eur
<https://github.com/PyPSA/pypsa-eur>`_, which clusters down the
electricity transmission substations in each European country based on
the k-means algorithm. This gives nodes which correspond to major load
and generation centres (typically cities).

The total number of nodes for Europe is set in the ``config.yaml`` file
under ``clusters``. The number of nodes can vary between 37, the number
of independent countries / synchronous areas, and several
hundred. With 200-300 nodes the model needs 100-150 GB RAM to solve
with a commerical solver like Gurobi.


Not all of the sectors are at the full nodal resolution, and some
demand for some sectors is distributed to nodes using heuristics that
need to be corrected. Some networks are copper-plated to reduce
computational times.

For example:

Electricity network: nodal.

Electricity residential and commercial demand: nodal, distributed in
each country based on population and GDP.

Electricity demand in industry: based on the location of industrial
facilities from `HotMaps database <https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database>`_.

Building heating demand: nodal, distributed in each country based on
population.

Industry demand: nodal, distributed in each country based on
locations of industry from `HotMaps database <https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database>`_.

Hydrogen network: nodal.

Methane network: single node for Europe, since future demand is so
low and no bottlenecks are expected.

Solid biomass: choice between single node for Europe and nodal where biomass
potential is regionally disaggregated (currently given per country,
then distributed by population density within)
and transport of solid biomass is possible.

CO2:  single node for Europe, but a transport and storage cost is added for
sequestered CO2. Optionally: nodal, with CO2 transport via pipelines.

Liquid hydrocarbons: single node for Europe, since transport costs for
liquids are low.
