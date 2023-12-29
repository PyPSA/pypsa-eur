..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Validation
##########################################

The PyPSA-Eur model workflow provides a built-in mechanism for validation. This allows users to contrast the outcomes of network optimization against the historical behaviour of the European power system. The snakemake rule ``validate_elec_networks`` enables this by generating comparative figures that encapsulate key data points such as dispatch carrier, cross-border flows, and market prices per price zone.

These comparisons utilize data from the 2019 ENTSO-E Transparency Platform. To enable this, an ENTSO-E API key must be inserted into the ``config.yaml`` file. Detailed steps for this process can be found in the user guide `here <https://transparency.entsoe.eu/content/static_content/Static%20content/web%20api/Guide.html>`_.

Once the API key is set, the validation workflow can be triggered by running the following command:

    snakemake validate_elec_networks --configfile config/config.validation.yaml -c8


The configuration file `config/config.validation.yaml` contains the following parameters:

.. literalinclude:: ../config/config.validation.yaml
   :language: yaml

The setup uses monthly varying fuel prices for gas, lignite, coal and oil as well as CO2 prices, which are created by the script ``build_monthly_prices``. Upon completion of the validation process, the resulting network and generated figures will be stored in the ``results/validation`` directory for further analysis.


Results
=======

By the time of writing the comparison with the historical data shows partially accurate, partially improvable results. The following figures show the comparison of the dispatch of the different carriers.

.. image:: ../graphics/validation_seasonal_operation_area_elec_s_37_ec_lv1.0_Ept.png
   :width: 100%
   :align: center

.. image:: ../graphics/validation_production_bar_elec_s_37_ec_lv1.0_Ept.png
   :width: 100%
   :align: center



Issues and possible improvements
--------------------------------

**Overestimated dispatch of wind and solar:** Renewable potentials of wind and solar are slightly overestimated in the model. This leads to a higher dispatch of these carriers than in the historical data. In particular, the solar dispatch during winter is overestimated.

**Coal - Lignite fuel switch:** The model has a fuel switch from coal to lignite. This might result from non-captured subsidies for lignite and coal in the model. In order to fix the fuel switch from coal to lignite, a manual cost correction was added to the script ``build_monthly_prices``.

**Planned outages of nuclear power plants:** Planned outages of nuclear power plants are not captured in the model. This leads to a underestimated dispatch of nuclear power plants in winter and a overestimated dispatch in summer. This point is hard to fix, since the planned outages are not published in the ENTSO-E Transparency Platform.

**False classification of run-of-river power plants:** Some run-of-river power plants are classified as hydro power plants in the model. This leads to a general overestimation of the hydro power dispatch. In particular, Swedish hydro power plants are overestimated.

**Load shedding:** Due to constraint NTC's (crossborder capacities), the model has to shed load in some regions. This leads to a high market prices in the regions which drive the average market price up. Further fine-tuning of the NTC's is needed to avoid load shedding.
