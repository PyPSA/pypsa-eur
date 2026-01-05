# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Industry configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#industry
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class IndustryConfig(ConfigModel):
    """Configuration for `industry` settings."""

    St_primary_fraction: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.6,
            2025: 0.55,
            2030: 0.5,
            2035: 0.45,
            2040: 0.4,
            2045: 0.35,
            2050: 0.3,
        },
        description="The fraction of steel produced via primary route versus secondary route (scrap+EAF). Current fraction is 0.6.",
    )
    DRI_fraction: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0,
            2025: 0,
            2030: 0.05,
            2035: 0.2,
            2040: 0.4,
            2045: 0.7,
            2050: 1,
        },
        description="The fraction of the primary route DRI + EAF.",
    )
    H2_DRI: float = Field(
        1.7,
        description="The hydrogen consumption in Direct Reduced Iron (DRI) Mwh_H2 LHV/ton_Steel from 51kgH2/tSt in `Vogl et al (2018) <https://doi.org/10.1016/j.jclepro.2018.08.279>`_.",
    )
    elec_DRI: float = Field(
        0.322,
        description="The electricity consumed in Direct Reduced Iron (DRI) shaft. From `HYBRIT brochure <https://ssabwebsitecdn.azureedge.net/-/media/hybrit/files/hybrit_brochure.pdf>`_.",
    )
    Al_primary_fraction: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.4,
            2025: 0.375,
            2030: 0.35,
            2035: 0.325,
            2040: 0.3,
            2045: 0.25,
            2050: 0.2,
        },
        description="The fraction of aluminium produced via the primary route versus scrap. Current fraction is 0.4.",
    )
    MWh_NH3_per_tNH3: float = Field(
        5.166,
        description="The energy amount per ton of ammonia (LHV).",
    )
    MWh_CH4_per_tNH3_SMR: float = Field(
        10.8,
        description="The energy amount of methane needed to produce a ton of ammonia using steam methane reforming (SMR). Value derived from 2012's demand from `Center for European Policy Studies (2008) <https://ec.europa.eu/docsroom/documents/4165/attachments/1/translations/en/renditions/pdf>`_.",
    )
    MWh_elec_per_tNH3_SMR: float = Field(
        0.7,
        description="The energy amount of electricity needed to produce a ton of ammonia using steam methane reforming (SMR). same source, assuming 94-6% split methane-elec of total energy demand 11.5 MWh/tNH3.",
    )
    MWh_H2_per_tNH3_electrolysis: float = Field(
        5.93,
        description="The energy amount of hydrogen needed to produce a ton of ammonia using Haber–Bosch process. From `Wang et al (2018) <https://doi.org/10.1016/j.joule.2018.04.017>`_, Base value assumed around 0.197 tH2/tHN3 (>3/17 since some H2 lost and used for energy).",
    )
    MWh_elec_per_tNH3_electrolysis: float = Field(
        0.2473,
        description="The energy amount of electricity needed to produce a ton of ammonia using Haber–Bosch process. From `Wang et al (2018) <https://doi.org/10.1016/j.joule.2018.04.017>`_, Table 13 (air separation and HB).",
    )
    MWh_NH3_per_MWh_H2_cracker: float = Field(
        1.46,
        description="The energy amount of amonia needed to produce an energy amount hydrogen using ammonia cracker.",
    )
    NH3_process_emissions: float = Field(
        24.5,
        description="The emission of ammonia production from steam methane reforming (SMR). From UNFCCC for 2015 for EU28.",
    )
    petrochemical_process_emissions: float = Field(
        25.5,
        description="The emission of petrochemical production. From UNFCCC for 2015 for EU28.",
    )
    HVC_primary_fraction: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.88,
            2025: 0.85,
            2030: 0.78,
            2035: 0.7,
            2040: 0.6,
            2045: 0.5,
            2050: 0.4,
        },
        description="The fraction of high value chemicals (HVC) produced via primary route.",
    )
    HVC_mechanical_recycling_fraction: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.12,
            2025: 0.15,
            2030: 0.18,
            2035: 0.21,
            2040: 0.24,
            2045: 0.27,
            2050: 0.30,
        },
        description="The fraction of high value chemicals (HVC) produced using mechanical recycling.",
    )
    HVC_chemical_recycling_fraction: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.0,
            2025: 0.0,
            2030: 0.04,
            2035: 0.08,
            2040: 0.12,
            2045: 0.16,
            2050: 0.20,
        },
        description="The fraction of high value chemicals (HVC) produced using chemical recycling.",
    )
    HVC_environment_sequestration_fraction: float = Field(
        0.0,
        description="The fraction of high value chemicals (HVC) put into landfill resulting in additional carbon sequestration. The default value is 0.",
    )
    waste_to_energy: bool = Field(
        False,
        description="Switch to enable expansion of waste to energy CHPs for conversion of plastics. Default is false.",
    )
    waste_to_energy_cc: bool = Field(
        False,
        description="Switch to enable expansion of waste to energy CHPs for conversion of plastics with carbon capture. Default is false.",
    )
    sector_ratios_fraction_future: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.0,
            2025: 0.05,
            2030: 0.2,
            2035: 0.45,
            2040: 0.7,
            2045: 0.85,
            2050: 1.0,
        },
        description="The fraction of total progress in fuel and process switching achieved in the industry sector.",
    )
    basic_chemicals_without_NH3_production_today: float = Field(
        69.0,
        description="The amount of basic chemicals produced without ammonia (= 86 Mtethylene-equiv - 17 MtNH3).",
    )
    HVC_production_today: float = Field(
        52.0,
        description="The amount of high value chemicals (HVC) produced. This includes ethylene, propylene and BTX. From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, Figure 16, page 107.",
    )
    MWh_elec_per_tHVC_mechanical_recycling: float = Field(
        0.547,
        description="The energy amount of electricity needed to produce a ton of high value chemical (HVC) using mechanical recycling. From SI of `Meys et al (2020) <https://doi.org/10.1016/j.resconrec.2020.105010>`_, Table S5, for HDPE, PP, PS, PET. LDPE would be 0.756.",
    )
    MWh_elec_per_tHVC_chemical_recycling: float = Field(
        6.9,
        description="The energy amount of electricity needed to produce a ton of high value chemical (HVC) using chemical recycling. The default value is based on pyrolysis and electric steam cracking. From `Material Economics (2019) <https://materialeconomics.com/latest-updates/industrial-transformation-2050>`_, page 125.",
    )
    chlorine_production_today: float = Field(
        9.58,
        description="The amount of chlorine produced. From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, Table 7, page 43.",
    )
    MWh_elec_per_tCl: float = Field(
        3.6,
        description="The energy amount of electricity needed to produce a ton of chlorine. From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, Table 6 page 43.",
    )
    MWh_H2_per_tCl: float = Field(
        -0.9372,
        description="The energy amount of hydrogen needed to produce a ton of chlorine. The value is negative since hydrogen produced in chloralkali process. From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, page 43.",
    )
    methanol_production_today: float = Field(
        1.5,
        description="The amount of methanol produced. From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, page 62.",
    )
    MWh_elec_per_tMeOH: float = Field(
        0.167,
        description="The energy amount of electricity needed to produce a ton of methanol. From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, Table 14, page 65.",
    )
    MWh_CH4_per_tMeOH: float = Field(
        10.25,
        description="The energy amount of methane needed to produce a ton of methanol. From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, Table 14, page 65.",
    )
    MWh_MeOH_per_tMeOH: float = Field(
        5.528,
        description="The energy amount per ton of methanol (LHV). From `DECHEMA (2017) <https://dechema.de/dechema_media/Downloads/Positionspapiere/Technology_study_Low_carbon_energy_and_feedstock_for_the_European_chemical_industry-p-20002750.pdf>`_, page 74.",
    )
    hotmaps_locate_missing: bool = Field(
        False,
        description="Locate industrial sites without valid locations based on city and countries.",
    )
    reference_year: int = Field(
        2019,
        description="The year used as the baseline for industrial energy demand and production. Data extracted from `JRC-IDEES 2015 <https://data.jrc.ec.europa.eu/dataset/jrc-10110-10001>`_.",
    )
    oil_refining_emissions: float = Field(
        0.013,
        description="The emissions from oil fuel processing (e.g. oil in petrochemical refinieries). The default value of 0.013 tCO2/MWh is based on DE statistics for 2019; the EU value is very similar.",
    )
