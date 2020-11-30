.. _overnight:

##########################################
Overnight (greenfield) scenarios
##########################################

The default is to calculate a rebuilding of the energy system to meet demand, a so-called overnight or greenfield approach.

For this, use ``foresight : 'overnight'`` in ``config.yaml``, like the example in ``config.default.yaml``.

In this case, the ``planning_horizons : [2030]`` scenario parameter can be set to use the year from which cost and other technology assumptions are set (forecasts for 2030 in this case).
