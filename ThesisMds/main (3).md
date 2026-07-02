::: titlepage
[Master Thesis]{.smallcaps}\

![image](./Figures/dtu_logo.png)\

------------------------------------------------------------------------

\
Cost-Optimal Deployment of Carbon Dioxide Removal Technologies under
Standalone Carbon Removal Markets in Europe\

------------------------------------------------------------------------

\

  -------------------- ---------------
  **Name**             **Study nr.**
  Philipp Hock         s240459
  Felicia Theilacker   s240026
                       
                       
                       
  -------------------- ---------------

2026-07-02\
:::

# Abstract {#abstract .unnumbered}

Can a market-based carbon dioxide removal (CDR) credit drive enough
BECCS and DACCS deployment to meet European CDR targets? This study
answers that question by extending PyPSA-Eur with a standalone,
eligibility-aware CDR credit market, spatially differentiated geological
storage costs, competing DACCS variants, weather-dependent capture
energy demand, and low-, medium-, and high-cost technology pathways. The
model traces how BECCS and DACCS respond to credit prices between 2030
and 2050, identifies the target-binding credit price, and compares it
with buyer willingness-to-pay benchmarks from stated and revealed market
evidence.

The results show that CDR does not emerge at target scale under ETS-only
pricing. A dedicated credit of roughly 100--200 €/tCO$_2$ is required to
meet European CDR targets across years and cost sensitivities.
Deployment shifts from BECCS in 2030 toward DACCS dominance by 2040 and
2050, while capture and storage concentrate in a few storage-rich
countries, especially around the North Sea. Despite large additional
electricity and heat demand, average system-price effects remain modest
because new renewable generation and heat pumps supply most of the added
demand. On the demand side, favourable willingness-to-pay benchmarks can
cover the per-tonne credit price, but the wider annualised system cost
remains unevenly distributed across countries. A standalone CDR credit
can therefore mobilise deployment, but infrastructure coordination and
burden-sharing instruments are still needed to finance the system-level
transition.

# Acknowledgements {#acknowledgements .unnumbered}

We thank our supervisors, Claire-Marie Bergaentzlé and Mak Dukan, for
their guidance and close feedback throughout this work.

We also thank our friends and family for their support over the past
months.

# List of Key Abbreviations {#list-of-key-abbreviations .unnumbered}

::: acronym
:::

# Code and Data Availability {#code-and-data-availability .unnumbered}

The full model code, including the PyPSA-Eur extensions developed in
this study, is openly available at
<https://doi.org/10.5281/zenodo.20775742> [@theilacker_pypsa-eur_2026].
It builds on PyPSA-Eur [@horsch_pypsa-eur_2018], available under the MIT
License.

# Introduction

The European Union aims for climate neutrality by 2050 and a 55% net
reduction in greenhouse gas emissions by 2030 compared to 1990 levels
[@european_commission_european_2019]. Achieving these goals requires
economy-wide emission reductions, but some residual emissions will
persist. Sectors like industry, aviation, and agriculture cannot fully
decarbonize through electrification, efficiency, or fuel switching alone
[@morris_net_2023; @rosenow_europe_2022].

Carbon dioxide removal (CDR) complements, but does not replace, emission
reductions [@ipcc_carbon_2022]. This study examines Bioenergy with
Carbon Capture and Storage (BECCS) and Direct Air Carbon Capture and
Storage (DACCS), the CDR technologies most integrated with the energy
system and most frequently included in EU energy-climate scenarios
[@shukla_climate_2022; @minx_state_2024; @lux_potentials_2023; @abegg_expert_2024].

BECCS combines biomass with carbon capture and storage. Biomass absorbs
CO~2~ during growth. When it is converted for energy, the released
carbon is captured and stored underground. The result is net removal.
The availability of sustainable biomass limits BECCS's scale.

DACCS removes CO~2~ directly from ambient air using chemical sorbents,
with the captured CO~2~ compressed and stored underground. There are two
main DACCS processes: solid-sorbent (S-DAC), which regenerates at low
temperatures and can use waste or district heat, and liquid-sorbent
(L-DAC), which regenerates at high temperatures and primarily uses
electricity. Both approaches are energy intensive and currently limited
by cost and technological maturity
[@shukla_climate_2022; @farzan_kazemifar_review_2021].

In addition to BECCS and DACCS, carbon capture can be applied to fossil
and process emissions. This approach prevents new emissions rather than
removing existing atmospheric carbon, so it is considered abatement
rather than CDR. All capture methods compete for geological storage,
with costs and capacity varying by location.

Captured can also serve as feedstock for synthetic fuels. This thesis
classifies this fate as utilisation, not credited CDR, since the carbon
is re-emitted when the fuel is used. Only atmospheric or biogenic that
is geologically stored counts toward the CDR target.

The CDR volumes implied by 1.5°C pathways exceed what current plans
deliver, yet existing incentives target abatement rather than removal
[@minx_state_2024; @lamb_carbon_2024]. The EU Emissions Trading System
(EU ETS), the Carbon Border Adjustment Mechanism, and corporate net-zero
pledges all put a cost on emitting or reward emission cuts. None of them
pays for carbon taken back out of the atmosphere
[@rickels_integrating_2021; @johnstone_carbon_2026].

Demand for removal is instead emerging through the voluntary carbon
market (VCM). The VCM mirrors the ETS but trades the opposite unit: the
ETS sells allowances that each permit one tonne emitted, while the VCM
sells credits that each certify one tonne removed
[@dawes_voluntary_2023]. Companies buy these credits to offset emissions
they are not required to abate. Most trades are bilateral and private,
so the market stays small and lacks a common reference price
[@roston_carbon_2026]. Eligibility depends on carbon origin: the
biogenic share of BECCS and all DACCS qualify, while fossil and process
CO~2~ do not
[@european_parliament_and_council_of_the_european_union_regulation_2024].

Policy recognition of CDR is increasing. The EU Taxonomy and the Carbon
Removal Certification Framework provide a legal basis for removal but do
not assign a price. Proposed solutions such as reverse auctions,
contracts for difference, and ETS integration could establish a price
[@rickels_integrating_2021; @yang_policy_2024; @lundberg_missing_2022; @fridahl_potential_2024].
While these instruments vary in design, each ultimately sets a price per
tonne removed. This study focuses on that common element.
Figure [1](#fig:cdr_overview){reference-type="ref"
reference="fig:cdr_overview"} compares the two price signals: the
voluntary credit pays for removal, while the EU ETS allows only fossil
and process capture to avoid the emissions price. Both approaches
compete for the same geological storage.

<figure id="fig:cdr_overview" data-latex-placement="H">

<figcaption>Capture routes and price signals in this study. Only
biogenic or atmospheric CO<sub>2</sub> that is geologically stored
counts as credited CDR.</figcaption>
</figure>

Existing energy-system models usually deploy CDR to satisfy a climate
constraint. In those models, the removal volume is fixed by a carbon
budget, net-zero target, or removal target, and the reported carbon
price is the shadow value of meeting that constraint
[@xie_negative_2025]. Emerging CDR markets work the other way around:
developers receive a payment per tonne removed, and the uncertain
outcome is how much removal that payment induces market-wide. Existing
models therefore answer a planning question, but not the market question
addressed here. They also simplify where and at what cost removal
happens, often using uniform geological storage costs and representing
DACCS as a single technology with fixed energy demand
[@lux_potentials_2023; @okosun_global_2026; @fernandes_exploring_2026].
As a result, they do not fully resolve how a standalone credit price
would shape CDR volume, technology choice, location, and system cost
across Europe.

On the demand side, willingness-to-pay (WTP) estimates are based on
surveys and expert input, but no study compares WTP to a modelled,
cost-optimal supply cost
[@cdrfyi_pricing_2026; @yang_policy_2024; @oh_review_2025]. The broader
system costs of CDR, such as added competition for generation, grid,
heat, and biomass, are also not assessed in one spatially resolved
European model
[@carbon_gap_envisioning_2024; @ferrier_defossilizing_2025; @millinger_are_2022].

This study uses PyPSA-Eur to assess how CDR responds to a standalone
credit price in Europe. Standalone refers to a credit that operates as
an independent price mechanism for CDR, separate from the EU ETS. The
study examines how much CDR is deployed, which technologies are
selected, where capture and storage locate, and what system cost
follows.

The analysis uses two linked metrics. The first is the target-binding
credit price: the lowest tested credit price at which cost-optimal
deployment meets Europe's CDR target in a given year. The second is the
financing gap, measured in two parts. The price gap compares the
target-binding credit price with buyer willingness to pay (WTP), drawn
from survey and transaction evidence. The system-level gap measures the
additional annualised system cost of meeting the target relative to the
zero-credit baseline, including how that cost is distributed across
countries.

This study provides an initial assessment of the value of an EU-wide CDR
credit. It identifies the credit price per tonne needed to meet Europe's
CDR targets and the remaining cost that must be addressed by other
policy instruments.

## Research questions {#sec:researchquestions}

This leads to the following research question and sub-questions

*RQ: How do CDR credit prices shape the cost-optimal deployment of CDR
in Europe, and what financing gap remains?*

To answer this overarching question, the analysis is guided by the
following sub-questions (SQ):

- *SQ1: At what credit price does CDR deploy to meet Europe's CDR
  target, where do CDR and storage locate, and with what effect on the
  energy system?*

- *SQ2: Does buyer willingness to pay cover the target-binding price,
  and what system cost does the targeted CDR add across countries?*

## Contributions

- **Methodological.** CDR is modelled as a response to a standalone
  per-tonne credit price rather than to the shadow price of a binding
  cap. This study extends PyPSA-Eur with country-differentiated storage
  costs, two competing DAC variants, and weather-dependent capture
  energy, and repeats the analysis across a range of technology costs. A
  willingness-to-pay benchmark built from survey and transaction data is
  set against the resulting target-binding price, and the annualised
  system cost of reaching the target is measured separately.

- **Policy.** The model outputs are read against the existing CDR policy
  literature. The aim is to indicate what a standalone credit price can
  and cannot achieve, and where dedicated support might still be needed.
  This gives policymakers a reference point for how far a single price
  instrument can carry removal before other measures must take over.

## Structure of the thesis {#structure-of-the-thesis .unnumbered}

The remainder of this thesis is structured as follows. Section
[2](#sec:litreview){reference-type="ref" reference="sec:litreview"}
reviews the literature on CDR technology characteristics, energy system
modelling, and the CDR financing gap. Section
[3](#sec:methodology){reference-type="ref" reference="sec:methodology"}
outlines the methodological framework, including the supply-side model,
demand-side willingness-to-pay assessment, and gap analysis.
Section [4](#sec:results){reference-type="ref" reference="sec:results"}
presents modelling results across scenarios and years.
Section [5](#sec:discussion){reference-type="ref"
reference="sec:discussion"} discusses the findings in the context of
European CDR policy and market development.
Section [6](#sec:conclusion){reference-type="ref"
reference="sec:conclusion"} offers policy implications and directions
for future research.

# Literature review {#sec:litreview}

## Representation of CDR in energy system models {#sec:lr-cdr-esm}

Energy system models (ESMs) provide a structured framework for analysing
interactions among carbon removal technologies, energy supply,
infrastructure constraints, and price signals. Yet most European CDR
models deploy removal to meet climate targets. A market payment does not
enter the deployment decision. Modelling CDR as a market therefore
requires three elements that remain underdeveloped: price-responsive
deployment, spatially resolved technology and storage costs, and
robustness to technology uncertainty.

### Constraint-driven deployment versus credit-price deployment {#sec:lr-endogenous}

Most CDR modelling treats removal as the least-cost response to a
climate constraint, rather than as a market response to a credit price.
ESMs offer detailed technological and spatial modelling, while
integrated assessment models (IAMs) address the economy, energy, and
climate at a higher level of aggregation. In the studies reviewed by
@xie_negative_2025, removal volumes are determined by a carbon budget,
net-zero year, or temperature target. This approach is suitable for
analysing how to meet a fixed climate outcome, but less appropriate for
testing a CDR market, where the payment per tonne is fixed and the
deployed volume is the result. In constraint-driven models, removal is
required and the carbon price is the shadow value of that requirement.
In market-based models, removal is incentivised, and the key uncertainty
is whether the incentive is sufficient to drive deployment.

This constraint-driven approach also dominates the European studies most
relevant to this thesis. Six reviewed studies resolve Europe at the
country level
[@solano_rodriguez_decarbonizing_2016; @tatarewicz_role_2021; @seck_hydrogen_2022; @lux_potentials_2023; @rodrigues_narrative-driven_2022; @negri_tailored_2024],
but all are driven by binding targets rather than a standalone credit
price. @seck_hydrogen_2022 use the MIRET-EU model to simulate EU
carbon-neutrality pathways under the European Climate Law's binding 2050
net-zero target. BECCS and DACCS are deployed to offset residual
emissions, with deployment levels determined by the net-zero constraint,
the CO~2~ storage injection cap, and biomass availability, not by a
removal price. @negri_tailored_2024 introduce a stochastic-programming
method for electricity demand uncertainty and optimise BECCS and DACCS
across the 27 EU states and the UK to meet a mandated cumulative removal
target. @rodrigues_narrative-driven_2022 compare three ESMs, REMIND-EU,
PRIMES, and ETM-UCL, against a single European carbon budget for 2020 to
2050. Although the implied carbon price and technology mix vary across
models, CDR is still deployed as a least-cost response to the budget.

This pattern continues in the most recent and highest-resolution
European study. @fernandes_exploring_2026 extend the PyPSA-Eur model, a
90-node, 3-hourly framework, to include a portfolio of removal options,
but deploy them only to meet a net-zero constraint. Consequently, the
reported carbon price is again the dual value of that constraint.

A few recent studies move beyond this approach, but none models a
standalone credit acting directly on BECCS and DACCS choices in a
spatially resolved European system. This standalone structure reflects
how removal incentives are now implemented in practice.
@morrow_gcam-cdr_2023 extend the GCAM ESM model with a module that
allows CDR demand to respond to policy mechanisms such as independent
CDR targets, interregional trade, and growth-rate controls, but the
underlying model is an aggregated IAM with limited energy system detail.
@bistline_value_2025 provide the closest precedent: a detailed US
electricity system model shows that voluntary carbon markets can reduce
decarbonisation costs by up to 50%, by substituting for costly direct
mitigation. They therefore demonstrate price-based deployment, but not
in the European, spatially resolved BECCS-DACCS setting addressed here.

European CDR markets are developing in this direction. Voluntary reverse
auctions, such as Sweden's Klimatpremie, already pay developers a price
per tonne removed [@lundberg_missing_2022; @fridahl_potential_2024], and
proposed compliance instruments, including integration into the EU ETS,
follow a similar approach [@sultani_sequencing_2024]. Similar
compliance-linked approaches are also being considered in the UK and
Japan [@cdrfyi_pricing_2026]. In both voluntary and proposed compliance
markets, the instrument is a price signal to developers, not a
system-wide carbon budget. A standalone credit price therefore best
represents current market operations [@johnstone_carbon_2026].
Accordingly, this study introduces an exogenous, standalone CDR credit
price into a spatially resolved European model. The credit price is the
only removal-specific revenue, and BECCS and DACCS deployment respond
directly to it.

### Limited spatial and operational resolution of BECCS and DACCS {#sec:lr-spatial}

Existing European models determine CDR location based only on renewable
energy and biomass availability, treating capture and storage costs as
spatially uniform. This approach omits the key cost factor that
determines where removal can occur: the cost of geological storage.
@lux_potentials_2023 and @fernandes_exploring_2026 both apply a flat 10
€/tCO~2~ storage cost across all regions and resolve storage only
through a national potential cap, not differentiated costs. None of the
reviewed studies introduces spatially differentiated sequestration
costs, so all miss an important driver of removal location. Since
captured CO~2~ must be stored, a uniform cost removes the main economic
reason for removal to concentrate in certain regions. As a result, these
models site removal based solely on resource availability.

In addition to storage cost, two technical characteristics of DACCS are
also held uniform across regions and have not been tested as siting
drivers. First, the energy required for capture depends on ambient
temperature and humidity, which vary across European climates
[@okosun_global_2026]. However, existing energy-system models apply a
single capture energy demand to all regions
[@lux_potentials_2023; @fernandes_exploring_2026]. Second, DACCS is
treated as a single technology, though solid-sorbent DAC (S-DAC)
regenerates at low temperatures and liquid-sorbent DAC (L-DAC) requires
high temperatures. The cost-optimal variant depends on the local supply
of low- versus high-grade heat and the associated renewables
[@okosun_global_2026].

A third operational characteristic is how the energy for capture is
supplied. DAC requires two energy carriers: electricity and, for the
solid-sorbent variant, low-grade heat [@okosun_global_2026]. Most
energy-system studies that co-optimise removal with supply treat both as
fixed inputs, so neither carrier is sourced endogenously.
@prado_assessing_2023 draws DAC power from the grid as a fixed
electricity demand, @lux_potentials_2023 converts regeneration heat into
an electricity equivalent, and @bistline_impact_2021 supplies it from
natural gas. In all cases, the energy-supply mix is fixed before
optimisation and enters as an input. @lehtveer_beccs_2021 take a
different approach: DACCS draws power from the modelled generation mix,
and this endogenous sourcing lowers system cost relative to BECCS
despite a higher levelised cost of carbon. However, that model covers
only electricity, in a single copperplate region, and treats DAC as one
technology. It remains untested whether sourcing both electricity and
low-grade heat within the model, across countries with different storage
and energy costs, affects which supply technologies are built and at
what system cost.

This thesis tests these siting mechanisms together. It differentiates
geological storage costs by country for each capture technology,
calibrated to the storage-type cost ranges of the
@zero_emissions_platform_costs_2011 and @martinez_castilla_clean_2025.
DACCS is separated into competing solid- and liquid-sorbent variants,
with weather-dependent, country-level capture energy based on
@okosun_global_2026. DAC electricity and low-grade heat are represented
as carriers the model supplies from competing options. This approach
tests whether storage access, climate-dependent DAC performance, or
endogenous heat and electricity sourcing affect CDR location, technology
deployment, and system cost.

### Uncertainty in technology development {#sec:lr-cost-uncertainty}

BECCS and DACCS are emerging technologies with uncertain future costs,
energy demand, and efficiency. While costs are expected to decline as
capacity grows and technologies mature, the extent of these reductions
is unknown. The relative cost-effectiveness of each technology depends
on how these factors evolve together
[@abegg_expert_2024; @victor_impact_2024]. @xie_negative_2025 identify
this uncertainty as a primary driver of results.

No study tests a full cost range across BECCS and both DAC variants
together. @lux_potentials_2023 come closest, presenting three DACCS cost
cases, but only for one DAC variant in a single year.
@fernandes_exploring_2026 adjust costs once, from a 2030 to a 2050 set.
@rodrigues_narrative-driven_2022 swap entire models, preventing
isolation of individual technologies. @seck_hydrogen_2022 do not vary
CDR costs at all. None of these studies shows how deployment and
technology mix shift when costs move from optimistic to conservative
scenarios.

This thesis applies low-, medium-, and high-cost sensitivities for
BECCS, S-DAC, and L-DAC for 2030, 2040, and 2050, based on the Danish
Energy Agency catalogue [@danish_energy_agency_technology_2024]. These
sensitivities represent plausible technology-development pathways: lower
costs indicate faster improvement, while higher costs indicate slower
progress. Within each sensitivity, investment cost, energy demand, and
efficiency move together, so each case reflects a coherent technology
pathway across all three parameters. This approach tests whether the
credit price required to meet the CDR target remains robust across
plausible technology futures or if market feasibility depends on
optimistic cost assumptions.

## Market feasibility: price and system-cost gaps {#sec:lr-price-gap}

The modelling gaps above concern CDR supply: how much removal a credit
price induces, where it is built, and which technology is chosen. To
assess whether this supply is market-feasible, the target volume must be
linked to financing. Existing CDR-gap studies mainly establish the scale
of missing removal: @lamb_carbon_2024 compare IAM pathways with national
pledges and inventories, while @yang_global_2023 allocate required
removal across countries using land and storage capacity. These studies
show how much removal is missing, but not what it costs to close the gap
or how much of that cost a credit market could cover. This thesis does
not re-estimate the required removal volume. It takes the European
removal target as the quantity to be met and asks whether it can be
financed through the credit market. Overall, it focuses on two financing
questions: whether buyer WTP covers system-derived supply costs and what
residual system costs remain.

### The price gap: willingness to pay versus system-derived supply cost

The price-gap literature measures the per-tonne difference between
supplier costs and what a buyer will pay. It reads WTP two ways, each
with the opposite weakness.

The first approach, stated preference, uses survey evidence. The 2026
CDR.fyi and OPIS assessment finds 2030 supplier prices above buyer WTP
for both BECCS (\$279 versus \$217/tCO~2~) and DACCS (\$400 versus
\$361/tCO~2~) [@cdrfyi_pricing_2026]. @yang_policy_2024 provide similar
European survey evidence, but stated preferences cover a limited market
segment and may exceed realised prices [@rodemeier_willingness_2023].

The second approach, revealed preference, is based on prices buyers have
actually paid in voluntary market transactions. Although these
transactions are limited and concentrated among a few buyers, they
provide a market-based price signal alongside survey data
[@cdrfyi_pricing_2026; @johnstone_carbon_2026]. Policy benchmarks
support similar conclusions, whether from proposals to include removals
in the EU ETS [@rickels_integrating_2021] or from Sweden's reverse
auctions for BECCS [@lundberg_missing_2022; @fridahl_potential_2024]. In
all cases, supplier costs are determined externally through project
estimates, surveys, or expert judgment. Whether buyer payments can cover
supply costs derived from the energy system has not yet been tested.

Both stated- and revealed-preference evidence provide price benchmarks,
but neither offers a full price-quantity relationship for Europe. The
literature shows what some early buyers may pay per tonne for removal
but does not indicate how demand changes as prices rise. No study has
estimated a price-dependent European CDR demand curve. This thesis does
not attempt to model total demand. Instead, it tests whether the
model-derived credit price required to meet the CDR target falls within
stated and revealed willingness-to-pay benchmarks. The gap is measured
as a price gap, not as a complete demand-supply equilibrium.

### The system-cost gap: plant costs versus wider financing needs

The capital-gap literature estimates removal costs plant-by-plant rather
than system-wide. Policy and industry assessments multiply capacity
targets by per-tonne unit costs to produce a single European funding
estimate [@carbon_gap_envisioning_2024; @ferrier_defossilizing_2025].
Other studies use similar methods, including expert elicitation of BECCS
and DACCS cost trajectories [@abegg_expert_2024], interviews on expected
supplier revenues across countries [@kozian_global_2026], and reviews of
method-level marginal abatement costs
[@michaelowa_international_2023; @oh_review_2025].

These approaches miss the broader system impacts of scaling BECCS and
DACCS, including additional generation, heat, grid investment, and
biomass competition. @lehtveer_beccs_2021 show that CDR rankings can
change once flexibility and displaced generation are included, while
@millinger_are_2022 isolate biomass competition rather than full system
cost. @negri_tailored_2024 come closest, reporting the cost of a removal
target as a single European sum, but without country-level allocation or
comparison to a removal price. The annualised system cost of meeting a
CDR target, resolved by country and assessed separately from the
per-tonne credit price, has not been quantified.

This thesis reassesses the price and capital-gap literatures using costs
generated within the energy system model. For the price gap, it compares
the credit price required to meet the CDR target with stated and
revealed WTP benchmarks. For the system-cost gap, it calculates the
additional annualised system cost of achieving the target and treats
this as the residual financing need after credit revenues are
considered. The contribution is not another estimate of CDR potential,
but a test of market feasibility: whether credit prices can mobilise
supply, whether those prices align with buyer WTP, and what system costs
remain.

## Summary and research gap {#sec:lr-summary}

Five gaps in the literature motivate this thesis.

1.  **Shadow price, not a market signal.** Models set removal to hit a
    target, so their carbon price is its shadow price. A standalone
    credit stays untested.

2.  **Uniform storage and DACCS detail.** Storage cost is flat, DACCS is
    one variant, and its energy is fixed, so where removal locates is
    assumed, not solved.

3.  **Single technology pathway.** Studies run one cost pathway, so no
    result is shown to hold from optimistic to conservative costs.

4.  **Exogenous supply costs in the price gap.** WTP is set against
    survey and project costs, not a system where BECCS and DACCS compete
    for biomass, power, and storage.

5.  **Plant-level capital in the system-cost gap.** Capital is costed
    plant by plant, and the wider system cost is rarely estimated, never
    split by country or set against a price.

Taken together, these gaps leave the market feasibility of European CDR
untested. They map onto the two halves of the research question posed in
Section [1.1](#sec:researchquestions){reference-type="ref"
reference="sec:researchquestions"}. How credit prices shape deployment
requires a model where removal responds to a standalone credit (gap 1),
storage costs and DAC performance differ by country (gap 2), and results
hold across technology cost pathways (gap 3). What financing gap remains
requires comparing the required credit price with buyer WTP (gap 4) and
resolving the residual system cost by country (gap 5).

# Methodology {#sec:methodology}

## Modelling framework and key assumptions {#subsec:framework}

This study extends the open-source PyPSA-Eur model to test how CDR
credit prices drive CDR deployment in Europe. PyPSA-Eur resolves
electricity, storage, transmission, sector coupling, CO~2~ transport,
and geological storage across European nodes [@horsch_pypsa-eur_2018].
Captured can be routed between nodes through an extendable bidirectional
pipeline network, so capture and storage locations can differ. The
network formulation is described in
Appendix [7.1.10](#subsubsec:app_co2_transport){reference-type="ref"
reference="subsubsec:app_co2_transport"}. It finds the least-cost mix of
investment and dispatch under techno-economic, network, and policy
constraints [@neeraj_dhanraj_bokde_pypsa_2025]. The model runs at 96
nodes and uses myopic foresight, so 2030, 2040, and 2050 solve
sequentially, with earlier investments carried forward. The full
configuration is in
Appendix [7.1.1](#subsubsec:app_pypsa_config){reference-type="ref"
reference="subsubsec:app_pypsa_config"}; the model code is publicly
available [@theilacker_pypsa-eur_2026].

The methodology links CDR supply with market financing. PyPSA-Eur
estimates how much removal a given credit price delivers. External WTP
benchmarks then test whether buyers could cover the resulting
target-binding price.
Figure [2](#fig:method_schematic){reference-type="ref"
reference="fig:method_schematic"} shows how the supply-side outputs and
demand-side benchmarks meet in the financing-gap analysis.

<figure id="fig:method_schematic" data-latex-placement="H">

<figcaption>Methodology overview</figcaption>
</figure>

The supply-side handoff variable is the target-binding credit price,
$\mathrm{CP}^{*}$, the lowest tested credit price at which credited CDR
reaches the yearly target (defined formally in
Section [3.1.2](#subsubsec:supply-curve){reference-type="ref"
reference="subsubsec:supply-curve"}). Scenario 1 derives
$\mathrm{CP}^{*}$ from the modelled supply curve. Scenario 2 sets
$\mathrm{CP}^{*}$ against buyer willingness to pay (WTP) and sets the
target-binding system cost against the baseline.
Section [3.1.4](#subsubsec:gaps){reference-type="ref"
reference="subsubsec:gaps"} defines both gap measures. The core
assumptions used across the runs are summarised in Appendix
Table [7](#tab:core_assumptions){reference-type="ref"
reference="tab:core_assumptions"}.

### Extensions to PyPSA-Eur {#subsubsec:supply-modules}

The two model extensions developed for the study are detailed below:

#### Module 1: standalone CDR credit market.

The first module lets deployment respond to an explicit credit price.
Removal no longer depends only on a binding emissions cap. The credit
enters the objective as a revenue term for every eligible tonne durably
stored:

$$\begin{equation}
    C^{\mathrm{credit}} = -\,\mathrm{CP}\sum_{o \in \mathcal{E}} s^{\mathrm{cred}}_{o},
    \label{eq:credit}
\end{equation}$$ where $\mathrm{CP}$ is the credit price and
$s^{\mathrm{cred}}_{o}$ is the annual credited removal attributed to
eligible origin $o$ (per model year; the year index is suppressed for
readability), with $\mathcal{E}={\text{BECCS},\ \text{DACCS}}$.

$$\begin{equation}
    s^{\mathrm{cred}}_{o} \le K_{o}\quad \forall\, o \in \mathcal{E},
    \qquad
    \sum_{o \in \mathcal{E}} s^{\mathrm{cred}}_{o} \le S^{\mathrm{seq}},
    \qquad \mathcal{E} = \{\text{BECCS},\ \text{DACCS}\},
    \label{eq:eligible}
\end{equation}$$ where $K_{o}$ is the annual CO~2~ captured by eligible
origin $o$ and $S^{\mathrm{seq}}$ is the total physical geological
sequestration in that year. Fossil & process capture are excluded from
$\mathcal{E}$ and earn no credit.

The credit applies only to eligible removals, and only up to the amount
of CO~2~ captured and delivered to storage. Since captured CO~2~ is
combined in a single storage system, the allocation of
$S^{\mathrm{seq}}$ across eligible origins is not determined by the
optimization; it is assigned ex post using the attribution rule in
Appendix [7.1.4](#subsubsec:app_waterfall){reference-type="ref"
reference="subsubsec:app_waterfall"}. Accordingly,
$s^{\mathrm{cred}}_{o}$ is a bounded attributed quantity, not a
physically traced per-tonne flow. The credit market module is
implemented in `scripts/solve_network.py` (function
`add_cdr_credit_accounting`) of the published repository
[@theilacker_pypsa-eur_2026].

The second constraint prevents double compensation from the model's
existing carbon price. In PyPSA-Eur, the global `co2 atmosphere` bus has
marginal value $p^{\mathrm{ETS}}{t}$, so each tonne removed by DACCS or
BECCS is already rewarded at the prevailing ETS price. Without the
constraint, an eligible CDR tonne would receive both the carbon value
and the credit, i.e. $\mathrm{CP}+p^{\mathrm{ETS}}{t}$ per tonne.

To leave the credit as the only CDR-specific revenue, an offsetting cost
equal to that implicit value is charged on each eligible capture link:
$$\begin{equation}
    C^{\mathrm{cancel}} = \sum_{n,t} p^{\mathrm{ETS}}_{t}\, a_{n,t},
    \label{eq:cancel}
\end{equation}$$ with $a_{n,t}$ the atmospheric CO~2~ removed at node
$n$ in period $t$. The implicit reward and the added cost then cancel
term by term, $$\begin{equation}
    \underbrace{-\,p^{\mathrm{ETS}}_{t}\, a_{n,t}}_{\text{implicit ETS reward}}
    \;+\;
    \underbrace{p^{\mathrm{ETS}}_{t}\, a_{n,t}}_{\text{Eq.~\ref{eq:cancel}}}
    \;=\; 0,
    \label{eq:net_removal}
\end{equation}$$ so the only CDR-specific term left in the objective is
the credit of Eq. [\[eq:credit\]](#eq:credit){reference-type="ref"
reference="eq:credit"}. The cancellation term makes the CDR credit, not
the ETS price, the CDR-specific revenue. It is an accounting construct,
not a resource cost, so it is excluded wherever system cost is reported
(Eq. [\[eq:sys_cost\]](#eq:sys_cost){reference-type="ref"
reference="eq:sys_cost"}). The ETS cancellation term of
Eq. [\[eq:cancel\]](#eq:cancel){reference-type="ref"
reference="eq:cancel"} is applied at network construction
[@theilacker_pypsa-eur_2026], in (function ).

#### Module 2: spatially resolved CDR representation.

The second module resolves where CDR is built and which DAC technology
is selected at each node. Geological storage costs are differentiated by
country, $c^{\mathrm{store}}_{n}$, using formation-type data from the
Zero Emissions Platform and the JRC CETO database
(Appendix [7.1.9](#subsubsec:app_sequestration_costs_potential){reference-type="ref"
reference="subsubsec:app_sequestration_costs_potential"}).

The optimiser chooses between solid-sorbent (S-DAC) and liquid-sorbent
(L-DAC) at each node, given their energy-carrier needs, capital costs,
and efficiency paths
(Appendix [7.1.5](#subsubsec:app_daccs_variants){reference-type="ref"
reference="subsubsec:app_daccs_variants"}). DAC energy demand is scaled
by a country-level weather factor $\varphi_{n}$ from ERA5 temperature
and dew-point reanalysis, following @okosun_global_2026. S-DAC reject
heat is credited at 0.52 of its full district-heating value to reflect
its low temperature. The spatial CDR module is implemented in
`scripts/prepare_sector_network.py`: country-differentiated
sequestration costs in function `add_co2_tracking` and the two DAC
variants with weather-dependent energy demand in function `add_dac`
[@theilacker_pypsa-eur_2026].

### Supply-side credit-price sweep {#subsubsec:supply-curve}

For each model year and CDR cost sensitivity, the model is solved at
credit prices from 50 to 500, in 50 EUR/tCO~2~ steps. Each solve gives
the cost-optimal CDR volume at that credit price. These price-volume
points define the CDR supply curve for that year and cost sensitivity.

The target-binding credit price is the lowest tested credit price at
which credited CDR reaches the CDR target: $$\begin{equation}
\mathrm{CP}^{*}_{y,k}
=
\min \left\{
\mathrm{CP} :
R_{y,k}(\mathrm{CP}) \geq R^{\mathrm{target}}_{y}
\right\},
\label{eq:cp_star}
\end{equation}$$ where $y \in \{2030,2040,2050\}$ is the model year, $k$
is the CDR cost sensitivity, $R_{y,k}(\mathrm{CP})$ is credited CDR
deployment, and $R^{\mathrm{target}}_{y}$ is the CDR target. The adopted
targets are 22, 160, and 514 MtCO~2~/yr in 2030, 2040, and 2050, derived
from EU-31 and UK CDR removal targets as detailed in
Appendix [7.1.11](#subsubsec:app_policy_goals){reference-type="ref"
reference="subsubsec:app_policy_goals"}. The CDR targets are not imposed
as optimisation constraints. They are applied after the sweep to
identify $\mathrm{CP}^{*}$. The zero-credit baseline is solved
separately and is used as the reference for deployment and system-cost
deltas. The eleven-step credit-price configurations for each cost
sensitivity are generated from a single template by
`config/Thesis_Runs/generate_supply_curve_configs.py`, and the canonical
run configurations are archived under `config/Thesis_Runs/`
[@theilacker_pypsa-eur_2026].

### Demand-side: willingness to pay {#subsubsec:wtp}

Willingness to pay (WTP) is used as an external buyer-side benchmark for
credited CDR. It is compared with the target-binding credit price,
$\mathrm{CP}^{*}$, after the model runs. The benchmark is not a demand
curve: it tests whether buyer-side prices cover the required credit
price, not whether enough buyers exist at that price to purchase the
full target volume.

Two benchmarks are used. The stated benchmark comes from the 2026
CDR.fyi/OPIS buyer survey [@cdrfyi_pricing_2026]. The revealed benchmark
comes from executed Frontier and CDR.fyi transactions. For revealed
deals, the unit price is either disclosed directly or derived as
$$\begin{equation}
P_d = \frac{V_d}{Q_d},
\label{eq:deal_price}
\end{equation}$$ where $V_d$ is contract value and $Q_d$ is contracted
CDR volume. Deals without a public method and price basis are excluded.

For both benchmarks, BECCS and DACCS 2030 values are combined with equal
weight: $$\begin{equation}
\mathrm{WTP}_{2030}
=
\frac{1}{2}
\left(
\mathrm{WTP}^{\mathrm{BECCS}}_{2030}
+
\mathrm{WTP}^{\mathrm{DACCS}}_{2030}
\right).
\label{eq:wtp_central}
\end{equation}$$ Equal weighting keeps the demand benchmark independent
of the modelled technology mix and avoids volume-weighting DACCS out of
the early market evidence. All values are converted using 1 USD =
0.864 EUR.

The 2040 and 2050 values are projected from the 2030 anchors:
$$\begin{equation}
\mathrm{WTP}_{y}
=
\mathrm{WTP}_{2030}
(1+r)^{y-2030},
\label{eq:wtp_projection}
\end{equation}$$ where $r$ is the annual projection rate. These
later-year values are sensitivity assumptions, not forecasts.
Table [1](#tab:wtp_method_values){reference-type="ref"
reference="tab:wtp_method_values"} gives the benchmark values used in
the gap analysis.
Appendix [7.1.13](#app:wtp-construction){reference-type="ref"
reference="app:wtp-construction"} provides the detailed source values,
deal register, exclusion rules, and projection rationale.

::: {#tab:wtp_method_values}
  **Benchmark 2030**      **Low**   **Medium**   **High**   **Projection rates**
  ----------------------- --------- ------------ ---------- ----------------------------
  Stated survey                                             $-3.5\%,\ -1.0\%,\ +2.5\%$
  Revealed transactions                                     $-0.9\%,\ -1.1\%,\ +2.5\%$

  : WTP benchmark construction. Values are EUR/tCO~2~.
:::

### Measuring the gaps {#subsubsec:gaps}

The financing gap combines the supply-side target-binding credit price
from Section [3.1.2](#subsubsec:supply-curve){reference-type="ref"
reference="subsubsec:supply-curve"} with the WTP benchmark from
Section [3.1.3](#subsubsec:wtp){reference-type="ref"
reference="subsubsec:wtp"}. It is reported in two measures.

The price gap sets the supply-side price against the demand-side
benchmark: $$\begin{equation}
    \Delta G_{y} = \mathrm{CP}^{*}_{y} - \mathrm{WTP}_{y},
    \label{eq:price_gap}
\end{equation}$$ reported as a signed per-tonne value, negative where
WTP covers the price and positive where a gap remains.

The system-level gap stays on the supply side. It is the additional
annualised cost of reaching the target relative to the baseline
(Section [3.2](#subsec:scenarios){reference-type="ref"
reference="subsec:scenarios"}): $$\begin{equation}
    \Delta C_{y} = C^{\mathrm{corr}}_{\mathrm{target},y} - C^{\mathrm{corr}}_{\mathrm{base},y},
    \label{eq:sys_cost}
\end{equation}$$ where $C^{\mathrm{corr}}_{y}$ is the total annualised
system cost (capital plus real operating cost) with the ETS-cancellation
term (Eq. [\[eq:cancel\]](#eq:cancel){reference-type="ref"
reference="eq:cancel"}) removed, since no physical resource is consumed
when an allowance is retired. $\Delta C_{y}$ is decomposed into five
components (BECCS biomass, DACCS capture capital, CO~2~ transport and
storage, renewable generation, residual) and resolved at country level.

## Scenarios {#subsec:scenarios}

The analysis contains one baseline and two scenario steps, summarised in
Table [2](#tab:analysis_design){reference-type="ref"
reference="tab:analysis_design"}. Cost sensitivities affect the
supply-side model runs, while WTP sensitivities affect only the
post-processed price gap.

The baseline is the cost-optimal European system under ETS-only pricing,
with zero CDR credit price and medium technology costs. Scenario 1
traces the deployment response to the CDR credit-price sweep. Scenario 2
compares the resulting target-binding price and system cost with buyer
WTP and the baseline.

Cost sensitivity is applied to all model results through low-, medium-,
and high-technology-cost assumptions
[@danish_energy_agency_technology_2024]
(Appendix [7.1.6](#subsubsec:app_tech_cost_scenarios){reference-type="ref"
reference="subsubsec:app_tech_cost_scenarios"}). The three sensitivities
represent optimistic, central, and conservative technology-learning
outcomes, expressed as a lower or higher cost level. Each scales the
whole cost trajectory by a constant factor, so the relative fall over
time is identical across sensitivities and only the absolute level
changes. The three sensitivities represent optimistic, central, and
conservative technology-learning outcomes, expressed as a lower or
higher cost level. DAC investment costs fall by about 40% from 2030 to
2050, compared with 10% for BECCS and 16% for fossil & process capture.
DAC is therefore the main source of assumed learning, while the
point-source capture routes are treated as relatively mature.

The EU ETS price shapes the wider system, mainly because fossil &
process capture use the same geological storage as CDR. It does not
directly incentivise CDR removal because Module 1 cancels the ETS value
by construction. The central trajectory is 119, 279, and 463 EUR/tCO~2~
in 2030, 2040, and 2050 ([@gunther_carbon_2025], see Details in
Appendix [7.1.12](#subsubsec:app_ets_prices){reference-type="ref"
reference="subsubsec:app_ets_prices"}).
Appendix [7.2.4](#app:sensitivity_ets){reference-type="ref"
reference="app:sensitivity_ets"} tests low and high ETS paths.

::: {#tab:analysis_design}
+--------------+----------------------+----------+----------+-------------+----------------------------------+
|              |                      |          |          |             |                                  |
+:=============+:=====================+:=========+:=========+:============+:=================================+
|              | **Credit price**     | **CDR    | **WTP**  |             |                                  |
|              |                      | cost**   |          |             |                                  |
+--------------+----------------------+----------+----------+-------------+----------------------------------+
| **Baseline** | EUR/t                | medium   | not      | Zero-credit | - Removal technology mix         |
|              |                      |          | modelled | reference   |                                  |
|              |                      |          |          | system,     | -                                |
|              |                      |          |          | basis for   |                                  |
|              |                      |          |          | all         | - Spatial pattern                |
|              |                      |          |          | comparisons |                                  |
|              |                      |          |          |             | - Electricity & heat price       |
+--------------+----------------------+----------+----------+-------------+----------------------------------+
| **Scenario   | to 500 EUR/t, 50     | low,     | not      | Credit      | - [CP^\*^]{.underline} per cost  |
| 1:**\        | EUR/t steps          | medium,  | modelled | price swept |   sensitivity (incl. technology  |
| deployment   |                      | high     |          | until       |   mix)                           |
| response     |                      |          |          | credited    |                                  |
|              |                      |          |          | CDR meets   | -                                |
|              |                      |          |          | the CDR     |                                  |
|              |                      |          |          | target. The | - Spatial pattern                |
|              |                      |          |          | lowest such |                                  |
|              |                      |          |          | price is    | - Energy-system impact           |
|              |                      |          |          | CP\*.       |                                  |
|              |                      |          |          | Output      |                                  |
|              |                      |          |          | reported at |                                  |
|              |                      |          |          | CP\*        |                                  |
+--------------+----------------------+----------+----------+-------------+----------------------------------+
| **Scenario   | target-binding       | low,     | low,     | CP\* set    | - Price gap $=$                  |
| 2:**\        | [CP^\*^]{.underline} | medium,  | medium,  | against the |   [CP^\*^]{.underline}$\,-\,$WTP |
| financing    |                      | high     | high     | WTP bands,  |                                  |
| gap          |                      |          |          | and the     | - System-level gap $=$ $-$ , per |
|              |                      |          |          | Scenario 1  |   country                        |
|              |                      |          |          | system cost |                                  |
|              |                      |          |          | against the |                                  |
|              |                      |          |          | baseline,   |                                  |
|              |                      |          |          | to size the |                                  |
|              |                      |          |          | financing   |                                  |
|              |                      |          |          | gap         |                                  |
+--------------+----------------------+----------+----------+-------------+----------------------------------+

: Scenario design. All runs are solved for 2030, 2040 and 2050.
:::

Quantities carried between scenarios are flagged by underline style

*Disclaimer:* In this study, \"scenario\" denotes an analysis step.
Scenario 1 sweeps the credit price across cost sensitivities and
re-solves the model at each point. Scenario 2 runs no new optimisation
in PyPSA-Eur. It combines the Scenario 1 outputs ([CP^\*^]{.underline}
and system cost) with external WTP benchmarks and the baseline to size
the financing gap.

# Results {#sec:results}

The results are organised around the two research sub-questions (SQs).
The baseline establishes the zero-credit reference system. It is the
comparison point for all scenario outcomes: deployment volume,
technology mix, spatial pattern, and energy-system impacts. Scenario 1
then shows how target-binding CDR credit prices influence these factors.
Scenario 2 compares target-binding CDR credit prices with WTP benchmarks
and decomposes the system-cost delta needed to achieve those volumes.
Across the results, -related terms are kept separate as defined in Table
[3](#tab:co2-categories){reference-type="ref"
reference="tab:co2-categories"}.

::: {#tab:co2-categories}
  Term             Definition
  ---------------- --------------------------------------------------------------------
  CO~2~ Captured   All captured, from BECCS, DACCS, fossil & process CCS
  CDR              captured by removal technologies: BECCS and DACCS
  Credited CDR     CDR that is geologically stored and earns the credit
  Stored           All sent to geological storage: credited CDR, fossil & process CCS
  Utilisation      Captured routed into the synthetic-fuel pathway

  : Definitions of the -related terminology used across the results.
:::

## Baseline without a CDR credit price

The baseline delivers no target-scale CDR. The ETS price drives fossil &
process CCS and renewable expansion, but not atmospheric removal.

### Key metrics of the baseline model

The zero-credit baseline shows whether CDR appears without a CDR credit
price. It also defines the reference point against which Scenarios 1 and
2 are measured. The baseline uses medium cost sensitivity and a CDR
credit price of 0 /t. The EU ETS is then the only carbon price signal.
It follows the medium trajectory used throughout.
Table [4](#tab:baseline_metrics){reference-type="ref"
reference="tab:baseline_metrics"} summarises the main baseline outcomes.
Detailed figures are reported in
Appendix [7.2.1](#subsec:app_results_baseline){reference-type="ref"
reference="subsec:app_results_baseline"}.

The only CDR that appears is a modest amount of DACCS in 2050, equal to
11.3 Mt /yr, which remains far below the CDR target. Stored CO~2~ is
dominated by EU ETS-driven fossil & process CCS, peaking in 2040 before
declining in 2050 as emissions fall.

The average electricity price and heat price peak in 2040 rather than
falling monotonically with system costs. The shared 2040 peak is a
transitional pressure point. Demand growth, electrification, and
infrastructure expansion raise system costs in 2040. Further renewable
and heat-pump build-out then lowers prices by 2050.

::: {#tab:baseline_metrics}
  Metric                                  2030   2040   2050
  -------------------------------------- ------ ------ ------
  Physical geological storage (Mt /yr)     47     73     24
  of which fossil & process CCS            47     73     13
  of which DACCS                          0.0    0.0    11.3
  of which biogenic/BECCS                 0.0    0.0    0.0
  Average electricity price (€/MWh)        67     96     64
  Average heat price (€/MWh)               40     87     54
  Total system cost (bn €/yr)             877    675    638

  : Baseline headline metrics (medium-cost sensitivity, CDR credit price
  0 /t)
:::

The baseline system decarbonises mainly through renewable expansion and
declining fossil generation. Detailed generation-mix changes are
reported in
Appendix [7.2.1](#subsec:app_results_baseline){reference-type="ref"
reference="subsec:app_results_baseline"},
Figure [13](#fig:baseline_generation_mix){reference-type="ref"
reference="fig:baseline_generation_mix"}.

### Spatial pattern of baseline deployment {#sec:baseline:spatial}

The baseline shows where is captured and routed before CDR credit prices
affect deployment. Figure [3](#fig:baseline_map){reference-type="ref"
reference="fig:baseline_map"} reports country-level zero-credit flows
for the years 2030, 2040, and 2050, with detailed values in
Appendix [7.2.1](#subsec:app_results_baseline){reference-type="ref"
reference="subsec:app_results_baseline"},
Tables [19](#tab:app_baseline_spatial_2030){reference-type="ref"
reference="tab:app_baseline_spatial_2030"}--[21](#tab:app_baseline_spatial_2050){reference-type="ref"
reference="tab:app_baseline_spatial_2050"}. Each country node splits
captured CO~2~ by source on the left and its modelled fate, storage or
utilisation, on the right. Orange arrows show cross-border transport,
scaled by volume and direction. These results provide the spatial
reference point for SQ1.

capture is driven almost entirely by fossil & process CCS, not by CDR.
In 2030, the UK (22 Mt /yr, 48%), Norway (11 Mt /yr), and the
Netherlands (9 Mt /yr) account for most of the 47 Mt /yr captured, all
of it stored and none moved across borders. In 2040, capture reaches
87 Mt /yr and spreads to Italy, Spain, and Germany alongside the North
Sea countries, while storage peaks at 73 Mt /yr. The remaining 14 Mt /yr
is routed to synthetic fuels rather than stored, and inter-country
pipeline trade stays negligible (about 1 Mt /yr). In 2050, the baseline
routes most captured to utilisation instead of storage. Capture reaches
127 Mt /yr, while storage declines to 24 Mt /yr, as fossil & process is
used as methanol feedstock. This utilisation is concentrated in Spain,
Italy, and the UK, where captured CO~2~ aligns with low-cost renewable
hydrogen. The remaining storage is smaller and mainly distributed across
the Netherlands, Italy, Germany, and Spain. Cross-border transport
remains limited at about 7 Mt /yr. The only baseline CDR is 11.3 Mt /yr
from DACCS, which is well below the CDR target.

<figure id="fig:baseline_map" data-latex-placement="H">
<img src="./Figures/Results/Baseline/Baseline_spatial.png" />
<figcaption>Baseline capture fate, and cross-border transport at zero
CDR credit price (medium-cost sensitivity), 2030-2050.</figcaption>
</figure>

Deployed DAC is entirely S-DAC. L-DAC is available but never built, as
its higher electricity demand outweighs its lower capital cost at
European power prices
(Appendix [7.1.5](#subsubsec:app_daccs_variants){reference-type="ref"
reference="subsubsec:app_daccs_variants"}). All later references to DAC
denote S-DAC.

The baseline shows the ETS drives fossil & process CCS, but not
target-scale CDR.

## Deployment response to CDR credit prices {#sec:scenario1}

Scenario 1 answers SQ1: how the target-binding credit price affects
deployment volume, technology mix, spatial pattern, and energy-system
outcomes, all as deltas from the baseline.

### Target-binding credit prices {#sec:scenario1:volume}

The deployment response is reported at the target-binding credit price,
defined as the lowest price at which credited CDR meets the CDR target,
across all three cost sensitivities
(Figure [4](#fig:binding_credit_tech_split){reference-type="ref"
reference="fig:binding_credit_tech_split"}).

The target-binding credit price varies by at most one 50 /t  step within
each year. It is 100 /t  in 2030 for low- and medium-cost, and 150 /t 
for high-cost. In 2040, low-cost remains at 100 /t , while medium- and
high-cost require 150 /t . In 2050, low- and medium-cost clear at
100 /t , while high-cost reaches 150 /t .

Figure [4](#fig:binding_credit_tech_split){reference-type="ref"
reference="fig:binding_credit_tech_split"} presents CDR by technology,
illustrating the deployment response at the target-binding credit price.
Deployed CDR volume exceeds the target due to the 50€/t step size. The
figure also shows the total of credited CDR and CDR that is utilised but
not credited. Binding prices rise weakly with cost sensitivity and stay
within one 50 €/t step of 100 /t  in every year. They are not monotonic
over time. The medium-cost price rises to 150 /t  in 2040 and returns to
100 /t  in 2050, while the low- and high-cost prices stay flat. The
technology mix shifts from BECCS to DACCS. In 2030, deployment is
entirely BECCS, while DACCS surpasses BECCS in the 2040 medium-cost
scenario. BECCS does not scale further because biomass competition caps
it at about 70 Mt /yr. In 2050, DACCS accounts for 85.5-94.4% of CDR
across all cost sensitivities. Detailed values are provided in Appendix
[7.2.2](#app:scenario_tables){reference-type="ref"
reference="app:scenario_tables"},
Table [22](#tab:app_cost_headline){reference-type="ref"
reference="tab:app_cost_headline"}.

<figure id="fig:binding_credit_tech_split" data-latex-placement="H">
<img src="./Figures/Results/Scenario 1/binding_credit_tech_split.png" />
<figcaption>CDR deployment at the target-binding credit price, by year
and cost sensitivity, split into BECCS (green) and DACCS (blue). First
bar is baseline; dotted line marks CDR target</figcaption>
</figure>

### Spatial pattern of CDR deployment {#sec:scenario1:spatial}

At the target-binding credit price, the spatial results show the full
CO~2~ capture, storage, and utilisation response across Europe
(Figure [5](#fig:cdr_map){reference-type="ref"
reference="fig:cdr_map"}). Detailed values are reported in
Appendix [7.2.2](#app:scenario_tables){reference-type="ref"
reference="app:scenario_tables"}, Tables
[23](#tab:app_s1_spatial_2030){reference-type="ref"
reference="tab:app_s1_spatial_2030"},
[24](#tab:app_s1_spatial_2040){reference-type="ref"
reference="tab:app_s1_spatial_2040"},
[25](#tab:app_s1_spatial_2050){reference-type="ref"
reference="tab:app_s1_spatial_2050"}. Denmark and the UK lead in both
capture and storage, with storage confined to the North Sea. Storage
access and low-cost electricity drive deployment, while cross-border
transport and local climate matter little. Unless otherwise specified,
all figures refer to the medium-cost sensitivity.

The concentration of CDR increases over time. In 2030, all CDR is
derived from BECCS and distributed among biomass-rich countries. In
2040, DAC becomes dominant and capture becomes more concentrated:
Denmark alone hosts 45%, Denmark and the UK together 65%, and the top
five countries 80% of total capture. In 2050, these two hubs account for
70% of capture and 80% of storage. As a result, only a few countries
meet the European CDR target. Cost variations have a minor effect. Under
the high-cost sensitivity, capture is slightly more distributed: the top
two countries' share of 2050 capture falls from 70% to 63%, and DAC's
share of total capture drops from 79% to 73% as BECCS increases.

Capture and storage do not always occur in the same country. Spain and
Italy capture more than they store, while the Netherlands stores more
than it captures. However, most stored remains at the capture site, and
cross-border transport is limited to 11 Mt /yr in 2050, about 2% of
storage. The remaining difference between capture and storage is mainly
due to utilisation. In the medium-cost sensitivity, utilisation reaches
99 Mt /yr in 2050, supplied by fossil & process CCS.

Storage access and electricity prices have a greater impact on DAC
deployment than local climate, even though the model accounts for
climate-driven energy demand at the country level
(Appendix [7.1.7](#app:dac-weather){reference-type="ref"
reference="app:dac-weather"}). Per-tonne energy costs vary from 0.82
times the reference value in Norway to 1.06 in Greece, but this
variation does not correlate with DAC deployment ($\rho = 0.02$ in
2050). Instead, deployment aligns with storage access ($\rho = 0.84$ in
2040, $0.53$ in 2050). For example, Norway has Europe's most favorable
DAC climate and North Sea storage but captures only 16Mt, while Spain,
despite a less favorable climate, captures twice as much due to
inexpensive solar power. The climate-related cost difference is about
8,€/t,, or roughly one-sixth of a credit-price step. A dedicated
robustness check confirms that weather does not significantly affect
siting or CDR cost (Appendix
[7.2.2.5](#app:sensitivity_weather){reference-type="ref"
reference="app:sensitivity_weather"}).

<figure id="fig:cdr_map" data-latex-placement="H">
<img src="./Figures/Results/Scenario 1/Medium_spatial.png" />
<figcaption>CDR deployment at target-binding credit price, medium-cost
sensitivity, 2030-2050.</figcaption>
</figure>

### Energy system implications: demand, prices and supply mix {#sec:scenario1:energysystem}

At the target-binding credit price, the deployment response adds a large
energy-demand block, but system-wide electricity and heat price effects
remain small.
Figure [6](#fig:energy_demand_binding){reference-type="ref"
reference="fig:energy_demand_binding"} reports five effects, all
measured as deltas from the zero-credit baseline: CDR energy demand (a),
electricity price (b), electricity supply mix (c), heat price (d) and
heat supply mix (e). Detailed numbers are reported in
Appendix [7.2.2](#app:scenario_tables){reference-type="ref"
reference="app:scenario_tables"}.

The demand mix shifts after 2030 (Figure
[6](#fig:energy_demand_binding){reference-type="ref"
reference="fig:energy_demand_binding"} a). In 2030, all CDR is BECCS,
needing 77 to 142TWh/yr of biomass. From 2040, DAC enters every
sensitivity and draws both electricity and low-grade heat, with heat the
larger block (about 1.3 to 1.9 times the electricity load in 2050) and
met almost entirely by the same air-source heat pumps that decarbonise
buildings (Figure [6](#fig:energy_demand_binding){reference-type="ref"
reference="fig:energy_demand_binding"} e). BECCS biomass remains
relevant throughout.

The energy price effect stays modest and, for heat, even reverses as DAC
demand peaks (Figure
 [6](#fig:energy_demand_binding){reference-type="ref"
reference="fig:energy_demand_binding"} b--e). Electricity prices rise by
at most 4.0 €/MWh in 2040 and between only 2.7 and 5.2 €/MWh in 2050,
despite higher DAC demand, while heat prices rise between 3.8 and
6.5 €/MWh in 2040 but then fall between 2.2 and 6.0 €/MWh in 2050, when
DAC heat demand is largest. The reason sits on the supply side. New wind
covers most of the added electricity load (+114 to 211TWh/yr in 2050).
Air-source heat pumps running on that same wind supply almost all the
added heat (+134 to 386TWh/yr in 2050). Both loads are met by added
capacity rather than by higher prices on scarce generation. Gas and all
other sources change little (by at most 8 and 24 TWh/yr respectively).

The system-average price masks larger country-level changes. In the
medium 2050 sensitivity, the average electricity price rises by just
5.2 €/MWh, but the hosts diverge sharply: the Danish price falls by
27 €/MWh on added wind, the British price is unchanged, and the Spanish
price rises by 12 €/MWh. Heat shows the same spread, with the average
falling by 2.2 €/MWh in 2050 and the British district-heat price by
4.2 €/MWh. The 2030 deltas carry no DAC, since all CDR that year is
BECCS (Appendix [7.2.2.6](#app:host_price_effects){reference-type="ref"
reference="app:host_price_effects"} Tables
[26](#tab:app_host_elec_price_delta){reference-type="ref"
reference="tab:app_host_elec_price_delta"},
[27](#tab:app_host_heat_price_delta){reference-type="ref"
reference="tab:app_host_heat_price_delta"}).

<figure id="fig:energy_demand_binding" data-latex-placement="H">
<img src="./Figures/Results/Scenario 1/Scenario1_Energy_System.png" />
<figcaption>System deltas from the baseline under target-binding CDR.
Bars show low (L), medium (M) and high (H) cost. (a) CDR energy demand
by BECCS biomass, DAC electricity and DAC heat. (b) electricity price.
(c) electricity mix. (d) heat price. (e) heat mix.</figcaption>
</figure>

## Financing gap and system costs of the CDR target {#sec:scenario2}

Scenario 2 answers SQ2 with two gap measures: a price gap,
$\mathrm{CP}^{*}-\mathrm{WTP}$, and a system-level gap, the additional
annualised system cost of reaching the CDR target relative to the
baseline.

### WTP benchmark {#sec:wtp-results}

The benchmarks are buyer-side price signals. They show buyers' WTP for
durable CDR, which is lower than supplier breakeven costs
(Appendix [7.1.13](#app:wtp-construction){reference-type="ref"
reference="app:wtp-construction"}). The price is treated as uniform
across off-taker sectors. Two benchmarks are used. The stated benchmark
comes from survey responses. The revealed benchmark comes from executed
transactions. Table [5](#tab:wtp-benchmarks){reference-type="ref"
reference="tab:wtp-benchmarks"} reports both across three sensitivities
and three planning years.

::: {#tab:wtp-benchmarks}
+-------------+-------------------------------------------+---------------------------------------------------------+
|             | Stated preference (survey)                | Revealed preference (transactions)                      |
+:============+=========:+=========:+=========:+=========:+===========:+===========:+===========:+=================:+
| 2-5 (l)6-9  | 2030     | 2040^\*^ | 2050^\*^ | Rate     | 2030       | 2040^\*^   | 2050^\*^   | Rate             |
| Sensitivity |          |          |          |          |            |            |            |                  |
+-------------+----------+----------+----------+----------+------------+------------+------------+------------------+
| High        | 302.4    | 387.1    | 495.1    | $+2.5\%$ | 612        | 783        | 1,002      | $+2.5\%$         |
+-------------+----------+----------+----------+----------+------------+------------+------------+------------------+
| Medium      | 249.7    | 225.5    | 203.9    | $-1.0\%$ | 381        | 341        | 305        | $\approx -1.1\%$ |
+-------------+----------+----------+----------+----------+------------+------------+------------+------------------+
| Low         | 197.0    | 138.2    | 96.8     | $-3.5\%$ | 136        | 124        | 113        | $\approx -0.9\%$ |
+-------------+----------+----------+----------+----------+------------+------------+------------+------------------+

: WTP benchmarks for durable CDR across sensitivities and planning years
(€/tCO~2~).
:::

^\*^2040 and 2050 values are directional sensitivity assumptions; no
comparable data exist beyond 2030.

The revealed benchmark exceeds the stated one in the medium and high
sensitivities. That gap is a first-of-a-kind premium on early deals,
which surveys do not capture. The low sensitivity is the exception:
committed deals sit below stated intentions, consistent with buyers
contracting in line with what they reported. The high benchmarks stay
above the contemporaneous EU ETS price in every year, while the central
and low benchmarks decline over time and the high one rises. The 2040
and 2050 values are directional demand-side assumptions, not
observations.

### Price gap {#sec:scenario2:operator_gap}

The price gap is calculated by $\mathrm{CP}^{*}-\mathrm{WTP}$. Negative
values mean the benchmark WTP exceeds the target-binding credit price.
Positive values indicate an uncovered per-tonne gap.

Under the stated benchmark
(Figure [7](#fig:tier1_gap_stated){reference-type="ref"
reference="fig:tier1_gap_stated"}), the 2030 gap is closed across every
WTP and cost sensitivity: even the lowest stated WTP clears the highest
target-binding price. It stays closed through 2040 except under low WTP
at medium and high cost (12 /t). By 2050, falling WTP open a gap: high
and medium WTP stay covered, but low WTP leaves an uncovered gap of 4 to
54 /t across all cost sensitivities (Appendix
[7.2.3](#app:financing_gap){reference-type="ref"
reference="app:financing_gap"} Table
[28](#tab:app_gap_stated){reference-type="ref"
reference="tab:app_gap_stated"}).

<figure id="fig:tier1_gap_stated" data-latex-placement="H">
<img src="./Figures/Results/Scenario 2/tier1_operator_gap_stated.png" />
<figcaption>Price gap under the stated-preference WTP
benchmark.</figcaption>
</figure>

The revealed benchmark
(Figure [8](#fig:tier1_gap_revealed){reference-type="ref"
reference="fig:tier1_gap_revealed"}) opens the low-WTP gap earlier,
already in 2030 at high cost (14 /t). In 2040 the low-WTP gap covers
both medium and high cost (26 /t), and by 2050 it narrows back to high
cost alone (37 /t), leaving the medium-cost sensitivity covered. Under
medium and high WTP it stays closed every year (Appendix
[7.2.3](#app:financing_gap){reference-type="ref"
reference="app:financing_gap"} Table
[29](#tab:app_gap_revealed){reference-type="ref"
reference="tab:app_gap_revealed"}).

<figure id="fig:tier1_gap_revealed" data-latex-placement="H">
<img
src="./Figures/Results/Scenario 2/tier1_operator_gap_revealed.png" />
<figcaption>Price gap under the revealed-preference WTP
benchmark.</figcaption>
</figure>

Across both, medium and high WTP are always covered. Under low WTP, the
stated benchmark gaps from 2040 and fully by 2050, while the revealed
benchmark gaps earlier but narrows to high cost only by 2050.

### System-level gap and country-level distribution {#sec:scenario2:system_gap}

The system-level gap is the annual increase in system cost required to
meet the CDR target relative to the baseline
(Section [3.1.4](#subsubsec:gaps){reference-type="ref"
reference="subsubsec:gaps"}). It reports the full CDR deployment cost
including BECCS capital, biomass feedstock, DAC capital, transport and
storage, renewable generation, and residual fossil and other changes.

The total cost delta rises from 0.3-3.3 bn €/yr in 2030 to
15.0-32.4 bn €/yr in 2040 and 32.5-63.2 bn €/yr in 2050 across cost
sensitivities. In the medium sensitivity, it reaches 2.0, 23.5, and
46.0 bn €/yr in total.

Costs shift from BECCS capital in 2030 to DAC capital and transport and
storage in later years
(Figure [9](#fig:tier2_waterfall){reference-type="ref"
reference="fig:tier2_waterfall"}). In 2030, target-binding CDR consists
entirely of BECCS, making BECCS capital the main cost. In 2040, DACCS
enters and DAC capital, CO~2~ transport and storage, and renewable
generation become the main cost. In 2050 in the medium cost sensitivity,
DAC capital is the largest component at 18.4 bn €/yr, followed by
transport and storage at 15.9 bn €/yr, renewable generation at
6.9 bn €/yr, while BECCS and biomass feedstock are minor.

<figure id="fig:tier2_waterfall" data-latex-placement="H">
<img
src="./Figures/Results/Scenario 2/tier2_waterfall_medium_v12.png" />
<figcaption>Corrected system-level cost decomposition at the
target-binding credit price in the medium-cost sensitivity. Bars show
annualised deltas relative to the baseline.</figcaption>
</figure>

By 2050, transport and storage costs stay within 14.4-16.5 bn,€/yr
because stored volumes are similar, while remaining differences mainly
reflect DAC and BECCS capital costs, especially in the high-cost
sensitivity.

The spatial distribution varies by component, not just by total CDR
volume. Detailed per-country breakdowns for each year are provided in
Appendix [7.2.3](#app:financing_gap){reference-type="ref"
reference="app:financing_gap"}, Tables
[30](#tab:app_country_cost_2030){reference-type="ref"
reference="tab:app_country_cost_2030"},
[31](#tab:app_country_cost_2040){reference-type="ref"
reference="tab:app_country_cost_2040"},
[32](#tab:app_country_cost_2050){reference-type="ref"
reference="tab:app_country_cost_2050"}.

In 2030, the per-country cost increase comes almost entirely from BECCS
capital and its transport and storage, led by the UK at 0.6 bn €/yr and
Sweden at 0.4 bn €/yr. In 2040, costs shift toward DAC and storage:
Denmark reaches 10.5 bn €/yr, mainly from DAC and transport and storage,
while the UK reaches 5.1 bn €/yr with a larger renewable-generation
component. In 2050, the UK and Denmark dominate the cost geography, at
roughly 22 and 12 bn €/yr. The UK increase combines DAC capital,
transport and storage, and a large renewable-generation block, while
Denmark's is almost entirely DAC capital and storage. Spain follows at
about 4 bn €/yr, mainly from renewable generation.

Cost increases follow CDR hosting for the main hubs
(Figure [10](#fig:tier2_fairshare){reference-type="ref"
reference="fig:tier2_fairshare"}). Denmark and the UK see the largest
increases because they host most of DAC and storage infrastructure.
Outside these hubs, costs reflect the system response: BECCS costs
appear in biomass-rich countries such as Sweden, Italy, and Spain, while
renewable costs are distributed across Germany, Ireland, France,
Ukraine, and Spain. This is why countries hosting little CDR can show
positive cost deltas. France, for example, incurs 0.9 bn €/yr in 2050,
mainly from renewable generation and biomass costs, while hosting only
about 6 Mt /yr of CDR and no storage. Poland and Ukraine follow the same
pattern. Country totals do not sum to the European total because oil,
gas, and biomass are priced in a shared European pool, leaving about
2.2 bn €/yr unallocated in 2050.
Appendix [7.2.3](#app:financing_gap){reference-type="ref"
reference="app:financing_gap"} reconciles the difference. Negative
segments are savings, not missing costs. Because the model optimises
Europe as one system, adding the credit shifts some renewable, biomass,
and fossil activity across countries. Countries doing less of an
activity than in the baseline record negative segments. Italy in 2050 is
the clearest case, where cheaper power displaces local renewable
generation; in Norway and France, fossil activity shifts toward CDR
hubs.

<figure id="fig:tier2_fairshare" data-latex-placement="H">
<img src="./Figures/Results/Scenario 2/Scenario2_tier2_countries.png" />
<figcaption>Country-level system-cost delta versus the zero-credit
baseline by technology bucket (medium-cost sensitivity); top eight
countries per planning year.</figcaption>
</figure>

# Discussion {#sec:discussion}

The discussion separates the two problems raised by the results: what a
standalone CDR credit can mobilise, and what it leaves unresolved in
financing and cost allocation.

## Deployment price, location, and system effects {#sec:disc-trigger}

This subsection links the binding credit price, deployment geography,
technology mix, and energy-system response. Together, these results show
what a standalone CDR credit can mobilise before the financing-gap
analysis turns to affordability and cost allocation.

### Without a dedicated credit, there is no CDR, and 100-150€/t meets the CDR goal {#sec:DiscussionSQ1.1}

In the baseline, the ETS supports fossil & process CCS, but not
target-scale CDR. Fossil & process capture avoids emissions costs,
whereas CDR requires a dedicated removal payment. Once the CDR credit is
introduced, credited CDR reaches the target. The lowest credit price
that does so is the target-binding price. In 2050, the target is met at
100-150 €/t across cost sensitivities
(Section [4.2.1](#sec:scenario1:volume){reference-type="ref"
reference="sec:scenario1:volume"}).

The target-binding price is robust within the tested assumptions. Cost
sensitivities have only a limited effect on the target-binding price.
However, all cost sensitivities include learning over time, so CDR costs
still fall and deployment becomes more competitive in later years.
Across low-, medium-, and high-cost sensitivities, the price changes by
at most one 50 €/t step, mainly affecting the technology mix rather than
target achievement. The ETS robustness check shows the same pattern: the
binding credit remains 100 €/t in 2030 and 2050 across ETS paths and
differs by only one step in 2040, likely reflecting the price grid
rather than a systematic ETS effect. This limited interaction follows
from the accounting design, where credited CDR has first claim on
geological storage
(Appendix [7.1.4](#subsubsec:app_waterfall){reference-type="ref"
reference="subsubsec:app_waterfall"}).

This credit price reflects a market price, not the shadow price of a
binding constraint. Other European studies deploy CDR to meet a carbon
budget or net-zero target and report the marginal abatement cost
[@xie_negative_2025]. @fernandes_exploring_2026 model this in PyPSA-Eur
and find 723 €/t under climate neutrality. This is not comparable to the
100-150 €/t range here, as their figure includes both removal and
abatement, whereas this study prices only credited CDR. The difference
follows from the policy design and would remain even under identical CDR
cost assumptions. A shadow price applies a single signal to every tonne,
while a standalone credit isolates the price of CDR, consistent with
current European practice, where CDRs are outside the EU ETS.

Policy research shows that CDR requires a dedicated demand signal.
@yang_policy_2024 find that a removal industry does not form without
government intervention. @johnstone_carbon_2026 note that carbon markets
reward CDR only when it is accounted for separately from emission
reductions. This study quantifies that argument. The target-binding
price is the modelled clearing price: the minimum credit at which the
least-cost system reaches the 514 Mt goal, resolved to the nearest
50 €/t step. Because it excludes risk premia, financing frictions, and
certification costs, the 100-150 €/t estimate is a lower bound on the
realised credit price.

Two implications follow. First, learning lowers CDR costs over time, but
falling costs alone do not create deployment. Learning is therefore a
precondition for CDR to become cost-competitive, but not sufficient to
bring it forward. A dedicated demand signal is required at every tested
cost sensitivity, so policy cannot wait for cheaper capture before
beginning deployment. Cost reductions influence the technology mix, but
deployment starts only when the credit is introduced.

Second, the standalone credit separates CDR revenue from the ETS. In
this framework, the binding credit price is therefore largely
independent of the ETS path, even at 2050 scale. This result depends on
the modelled accounting rule and on the ETS being exogenous: a fully
responsive ETS that adjusts to large-scale CDR remains outside the scope
of this study (Section [5.4](#sec:limitations){reference-type="ref"
reference="sec:limitations"}).

### CDR concentrates in a few storage-rich countries, and the hosts shift as the credit moves from BECCS to DACCS

CO~2~ Capture concentration shifts from an ETS-driven capture geography
to a credit-driven CDR geography. Without the credit, CO~2~ capture
relies mainly on fossil & process CCS, located where residual emissions
and storage access overlap. With the credit, deployment aligns with
conditions for credited CDR. In 2030, BECCS can meet the smaller target
by leveraging existing biomass supply chains and lower initial costs.
However, it does not scale much further because biomass is also needed
for industrial heat, transport biofuels, and methanol production. In the
model, this competition limits BECCS to about 70 Mt/yr.

Beyond this scale, further CDR relies on DAC, sited by storage access
and low-cost energy rather than biomass availability. Credited CDR
therefore clusters around North Sea storage, led by Denmark and the UK,
with Norway, Germany, the Netherlands, Ireland, and Sweden forming a
secondary tier. As DAC costs decline, deployment shifts from dispersed
BECCS toward these storage and power hubs. The minor weather effect
reinforces this interpretation. Siting follows storage access and energy
cost, and a favourable capture climate has almost no influence on where
DAC is built. Utilisation creates a separate capture geography that does
not count as credited CDR. At the target-binding price, synthetic fuels
use DAC-derived from southern Europe, especially Spain and Italy. Spain
captures about 50 Mt/yr but stores only 8 Mt/yr, with the rest used in
synthetic fuels. Thus, capture shifts south while credited removals
remain tied to North Sea storage.

This BECCS-before-DACCS sequence aligns with @prado_assessing_2023,
though the underlying mechanism differs. While they impose a biomass
ceiling, here the ceiling results from competition with heat, CHP, and
industry. This competition is not only a model artefact:
@millinger_diversity_2025 and @fernandes_exploring_2026 also find that
biomass is directed toward higher-value uses, even when total supply is
greater. The choice of DAC variant depends on energy prices. The model
selects only S-DAC because, at European power prices, L-DAC's higher
electricity demand outweighs its lower capital cost, consistent with
@okosun_global_2026. Lower electricity prices could make L-DAC more
competitive.

The storage outcome depends more on access than on cost.
@lux_potentials_2023 use flat sequestration costs and place DACCS in
Sweden, Norway, Finland, Iberia, and the Baltic States, while
@fernandes_exploring_2026 use uniform costs capped by national potential
and still find concentration in the UK, Denmark, Italy, and Portugal.
This study identifies the UK and Denmark as the leading storage hosts,
with Norway, Germany, the Netherlands, Ireland, and Sweden forming a
secondary tier. The overlap indicates that storage availability drives
the main concentration, while differences in storage costs mainly
influence which secondary hosts are selected. Here, that secondary tier
stays in the North Sea region, consistent with observed European storage
development, where operational and near-term projects such as Northern
Lights and Porthos cluster there [@martinez_castilla_clean_2025]. The
large southern capture in Spain and Italy sit outside this storage
ranking because they serve utilisation rather than storage.

This suggests that future CDR targets are constrained more by DAC
deployment and storage access than by biomass supply. Policy should
focus on ensuring affordable storage and coordinating DAC expansion in
likely host regions, particularly the UK and Denmark. The
country-specific storage costs used here are stylised, so the result is
a regional ranking
(Appendix [7.1.9](#subsubsec:app_sequestration_costs_potential){reference-type="ref"
reference="subsubsec:app_sequestration_costs_potential"}). Future
research could refine this allocation using more detailed geological
data and real-world siting constraints.

### DAC adds large energy demand, but the model meets it through new supply

DAC scaling changes the energy-system response to the CDR target.
Achieving the 2050 goal adds electricity, low-temperature heat, and
biomass demand, with DAC heat as the largest component
(Section [4.2.3](#sec:scenario1:energysystem){reference-type="ref"
reference="sec:scenario1:energysystem"}). Because the selected S-DAC
pathway relies on low-temperature heat, the model meets this demand
mainly through new wind and heat-pump capacity, while gas generation
changes little. Average prices stay mostly stable because the model
meets the added demand with new capacity.

Local prices show this adjustment more clearly than the European
average, and they move in both directions. In the northern storage
hosts, added supply lowers prices: Danish electricity falls by 27 €/MWh
on new wind, and the UK district-heat price falls by 4.2 €/MWh on new
heat pumps
(Appendix [7.2.2.6](#app:host_price_effects){reference-type="ref"
reference="app:host_price_effects"}). In the south, where DAC is tied to
e-fuel production rather than to storage, the added demand outpaces
local supply expansion and the Spanish electricity price rises by about
12 €/MWh. These movements do not mean that DAC lowers real-world energy
costs. They show where concentrated DAC demand is absorbed by new power
and heat capacity, and where it instead tightens local prices. The model
allows this expansion freely, whereas real systems would face
permitting, grid, and acceptance constraints.

This result extends studies that co-optimise CDR and energy supply.
@lehtveer_beccs_2021, @bistline_impact_2021, and @prado_assessing_2023
show that CDR changes the generation mix, but often predefine the energy
source for DAC or focus on added capacity. Here, the model chooses the
supply response endogenously. The contribution is to show where the
adjustment appears: large-scale DAC is met through new wind and
heat-pump capacity in host regions rather than through higher European
average prices.

The practical implication is that DAC deployment must be planned
together with its energy supply. The modelled price stability depends on
frictionless capacity expansion at assumed costs, while real renewable,
heat-pump, and grid expansion face permitting, network, and acceptance
constraints
[@manocha_reducing_2025; @international_energy_agency_electricity_2026].
Planning for DAC has to include electricity, heat, and grid capacity
alongside capture and storage.

## Price-gap coverage and remaining system cost {#sec:disc-cover}

Two findings answer SQ2. The per-tonne price gap closes in 2030 under
favourable WTP sensitivities, but this does not prove that voluntary
demand scales. The system-level gap is different: the wider
infrastructure cost is not automatically allocated by a per-tonne
credit.

### Buyer WTP covers the target-binding price under medium and high sensitivities, while low WTP opens a gap by 2050 {#sec:disc-pricelevel}

In 2030, the credit price meets both WTP benchmarks, with the single
exception of the low revealed benchmark at high cost, where a small gap
remains
(Section [4.3.2](#sec:scenario2:operator_gap){reference-type="ref"
reference="sec:scenario2:operator_gap"}). This is a best-case reading.
The cost side is a floor, based on least-cost optimisation that excludes
risk premia, delays, supply-chain issues, and host-country financing.
The WTP side is a ceiling. It is an early-mover price applied to the
full target volume. Each adjustment would widen the gap, and they do not
offset one another. The nearly closed 2030 gap rests on the most
favourable assumptions on both sides and is an upper bound on
affordability.

The two WTP benchmarks are based on limited evidence. The stated
benchmark reflects survey intentions, whereas the revealed benchmark is
based on a small number of signed deals at an early-mover premium. The
relationship differs by sensitivity, revealed prices exceed stated
values in the medium and high sensitivities because early transactions
include first-of-a-kind premia
(Table [5](#tab:wtp-benchmarks){reference-type="ref"
reference="tab:wtp-benchmarks"}) Neither benchmark determines the price
buyers would pay at scale.

The financing-gap literature reports a gap. The 2026 CDR.fyi and OPIS
survey puts it at 62 USD/t for BECCS and 39 USD/t for DACCS
[@cdrfyi_pricing_2026], while this study finds the 2030 gap closed. The
two numbers measure different things. CDR.fyi subtracts WTP from a
supplier's asking price, which carries a profit margin and a specific
cost of capital. This study subtracts WTP from the model's credit price,
which has no margin and applies a 7% WACC.

One more simplification matters here. WTP is held flat across sectors
and CDR targets (Section [4.3.1](#sec:wtp-results){reference-type="ref"
reference="sec:wtp-results"}), but comes from a small market. Scaling to
the target shifts demand down the curve, where later buyers pay less
than the first. A flat WTP therefore overstates what a target-scale
market would bear. Durable demand remains with a few buyers
[@newclimate_institute_companies_2025], and signed offtakes fall below
survey levels as contract length increases (Appendix
Table [17](#tab:beccs-duration){reference-type="ref"
reference="tab:beccs-duration"}, [@frontier_climate_frontier_2025]).

On the most favourable reading, today's buyers can cover the 2030 price.
This should not be read as evidence that voluntary demand can scale.
Even with a cost floor and an inflated WTP benchmark, a gap opens in
later years and under the pessimistic benchmark
(Section [4.3.2](#sec:scenario2:operator_gap){reference-type="ref"
reference="sec:scenario2:operator_gap"}). The later-year values rely on
assumed WTP paths with no data beyond 2030, so they are illustrative
potential outcomes (Section [5.4](#sec:limitations){reference-type="ref"
reference="sec:limitations"}).

### Targeted CDR adds capture-chain and energy-system costs across host and non-host countries {#sec:disc-systemgap}

The system-level gap shows the infrastructure cost behind the stable
average prices reported in SQ1. In the model, DAC electricity and heat
demand are absorbed through new wind and heat-pump capacity rather than
through higher European average prices
(Section [4.2.3](#sec:scenario1:energysystem){reference-type="ref"
reference="sec:scenario1:energysystem"}). The gap therefore measures
where the system has to expand to meet the CDR target.

A per-tonne credit covers the capture chain more directly than the wider
system around it. Most of the additional cost in the UK and Denmark is
DAC capital, transport, and storage, which are tied to the CDR activity
(Appendix [32](#tab:app_country_cost_2050){reference-type="ref"
reference="tab:app_country_cost_2050"}). What the credit does not
directly allocate is the supporting infrastructure: renewable
generation, heat supply, and biomass use built into the wider system
response. This residual is small in 2030, when BECCS meets the target
with limited new build, but grows as DAC becomes dominant in 2040 and
2050.

The residual cost is unevenly distributed across three roles. CDR hosts
such as the UK and Denmark carry the largest increases because they
build capture, transport, storage, and supporting energy infrastructure.
Some countries, such as France, incur costs despite hosting little CDR
because the integrated model builds resources where they are cheapest.
Others benefit from reallocation and record lower attributed costs than
in the baseline. Paying capture operators is only part of the financing
problem. The remaining part is allocating system costs and benefits
across hosts, non-hosting payers, and beneficiaries.

Most financing-gap studies estimate CDR affordability from plant-level
costs, supplier asking prices, or projected revenues
[@carbon_gap_envisioning_2024; @ferrier_defossilizing_2025; @kozian_global_2026].
This is useful for assessing whether a tonne can be financed, but it
misses the wider system response. As @lehtveer_beccs_2021 note,
levelised CDR costs can be misleading because they exclude system-wide
effects. Studies that include those effects usually report them as an
aggregate system cost. @negri_tailored_2024, for example, discuss an EU
and UK system cost that is fairly shared across member states, but do
not allocate it by country. This study shows what that aggregation
hides: once CDR is embedded in a spatially resolved energy system, part
of the cost falls outside the capture unit and is not automatically
reached by a per-tonne credit.

This country allocation is still model-dependent. PyPSA-Eur assumes a
highly integrated European system, so it can place costs where they are
cheapest. Real grids are less integrated, and permitting, grid
expansion, and local acceptance would limit how this adjustment occurs
[@klopecka_electricity_2025]. Where supply cannot expand, the same DAC
demand would more likely appear as higher local prices. The country cost
figures therefore show where the CDR target adds investment and
operating cost within the modelled system.

The practical implication is that large-scale CDR requires a
cost-allocation mechanism in addition to CDR credit payments. A credit
can finance the captured tonne, but it does not by itself allocate
renewable, heat, biomass, and cross-border system costs across hosts,
payers, and beneficiaries. This is especially important because much of
the modelled hosting occurs in the UK, outside any burden-sharing rule
the EU could pass alone.

## Implications: a three-phase deployment path {#sec:disc-finance}

The findings imply that CDR policy should be sequenced because the
binding constraint changes with scale. In 2030, the main constraints are
verified early demand and limited BECCS supply. By 2040, the challenge
shifts to coordinated DAC, storage, power, and heat expansion. By 2050,
CDR is concentrated in a few host countries, so the unresolved problem
is cross-border cost allocation. Deployment therefore involves three
successive challenges, and a single credit price does not address all of
them.

This sequencing is consistent with the climate-policy literature, where
instruments vary with cost, scalability, and maturity
[@meckling_policy_2017; @schenuit_carbon_2021; @winkler_exploring_2025].
Table [6](#tab:deployment_phases){reference-type="ref"
reference="tab:deployment_phases"} links each phase to the instrument
that fits the constraint identified in the model.

In 2030, bilateral offtake through the VCM can fit small, early volumes
because it funds verified tonnes before compliance demand exists. It
cannot be assumed to scale: VCM remain fragmented, face quality
concerns, and lack mature price discovery
[@dawes_voluntary_2023; @roston_carbon_2026]. Early policy should
therefore focus on credible buyers, verified removals, and measurement,
reporting, and verification (MRV) systems before volumes grow
[@burke_conditional_2024]. Two modelling caveats remain: buyer volume is
assumed rather than modelled, and BECCS supply may not scale because
biomass competes with other uses.

The 2040 phase is dominated by coordination. Once biogenic capture
saturates, DACCS must scale sharply, while transport, storage, power,
and heat assets expand in step. Public procurement and contracts for
difference fit this phase because the state can set deployment volumes
and absorb early price risk until compliance demand forms. Reverse
auctions can procure the lowest-cost tonnes directly, while contracts
for difference provide revenue certainty by covering the gap to a
reference carbon price
[@lundberg_missing_2022; @winkler_exploring_2025]. Both shift early risk
to the state, matching developer calls for state support before
long-term compliance demand emerges [@yang_policy_2024].

By 2050, CDR concentrates in a few storage-rich countries
(Section [5.2.2](#sec:disc-systemgap){reference-type="ref"
reference="sec:disc-systemgap"}). A per-tonne credit pays the CDR
provider, but it does not allocate the wider system costs between host
and non-host countries. An ETS-style instrument becomes relevant because
it creates mandatory demand and distributes costs through allowance
allocation. Since much of the hosting is in the UK, such an instrument
would need a cross-border design that goes beyond EU internal rules.
This matters because a uniform subsidy does not distribute costs
uniformly [@franks_optimal_2023]. Integration is only credible once
certification, liability, and MRV rules are in place
[@burke_conditional_2024; @sultani_sequencing_2024].

::: {#tab:deployment_phases}
  **Year**   **Phase**   **Primary challenge**                                    **Proposed instrument**
  ---------- ----------- -------------------------------------------------------- ----------------------------------------------
  2030       Pilot       BECCS dominates and competes for scarce biomass          VCM and bilateral offtake
  2040       Ramp        Scaling capture, storage and supporting infrastructure   Public procurement, contracts for difference
  2050       Mature      CDR concentrated in a few host countries                 ETS integration with burden-sharing

  : The three deployment phases, their primary challenges, and the
  demand instrument suggested for each.
:::

These instrument assignments are qualitative fits to the modelled
constraints. The model does not test them. The main implication is that
demand creation and cost allocation are separate policy problems. VCM
and procurement can support early deployment, but a mature compliance
instrument is needed once residual system costs and cross-border
burden-sharing dominate.

## Limitations and future work {#sec:limitations}

The limitations are organised by the results they qualify: credit price,
spatial deployment, system impacts, and financing gaps.

### Limitations affecting the target-binding credit price

These results come from a least-cost model. Target-binding credit prices
are especially sensitive to price resolution, target assumptions, ETS
representation and market frictions.

The credit-price sweep uses 50 €/t increments from 50 to 500 €/t. Each
reported binding price is the lowest tested step that meets the target.
The true threshold can sit anywhere in the 50 €/t interval below it. As
a result, binding prices and the resulting price gaps are upper bounds
within the tested grid and may be biased upward.

CDR targets are treated as fixed policy quantities. The study does not
test alternative target trajectories. Different target levels would
change the point on the supply curve where $\mathrm{CP}^{*}$ is
determined and could alter the technology mix, storage pressure, and
financing gap.

The EU ETS price is treated as an exogenous path, not as a cap-and-trade
market with banking, allowance supply, or a Market Stability Reserve.
Because the CDR credit cancels the ETS value of removals, the binding
price does not respond to the ETS level in the tested years. Credited
CDR take first claim on geological storage, so higher fossil & process
CCS does not reduce the credited share. This independence holds within a
fixed-price ETS. A fully responsive ETS that adjusts to large-scale
capture could introduce coupling and is left to future work.

The target-binding price is also a modelled clearing price, not a
realised market price. The model assumes perfect competition, uniform
financing, and no project-development or certification frictions.
Realised credit prices would likely exceed the modelled thresholds.
Technology costs are represented by low-, medium-, and high-cost
sensitivities, so conclusions about limited price sensitivity hold only
within that range.

### Limitations affecting spatial deployment and system impacts

The location results offer regional cost rankings and do not identify
specific project sites. Storage costs and potentials are assessed at the
country level and do not account for project-level factors such as
geological quality, permitting, public acceptance, liability, port
access, or detailed pipeline constraints. The observed concentration in
storage-rich countries is a modelled system tendency. Actual siting will
also depend on the project-level factors listed above.

The BECCS limit results from modelled biomass competition. Uncertainties
remain around biomass potentials, sustainability constraints, trade, and
alternative high-value uses. The transition from BECCS to DACCS depends
on these assumptions.

Computational constraints influence the timing and location of
deployment. The study applies myopic foresight, solving each period
sequentially with brownfield inheritance and without considering future
periods. This method may understate early investments that gain value at
higher CDR volumes and may overstate the abruptness of DAC scaling
between 2040 and 2050. Spatial and temporal resolution are limited to
ensure model tractability.

The energy-system effects should be viewed as indicators of model
pressure. PyPSA-Eur enables coordinated European capacity expansion
under least-cost assumptions, whereas real systems encounter permitting,
grid, supply chain, and local acceptance constraints. Where supply
cannot expand as easily, concentrated DAC demand may lead to higher
local prices or delayed deployment. The model uses a single
representative weather year, which excludes year-to-year variability and
understates the seasonality of DAC capture energy and heat demand.

Carbon utilisation is included in the model but is only partially
activated. Methanol dominates among utilisation routes due to its
exogenous demand target, while other options remain at or near zero
without similar demand or utilisation-specific revenue. Since the CDR
credit applies only to eligible that is geologically stored, utilised
carbon does not earn CDR revenue and remains subject to the ETS price if
later re-emitted. The utilisation-versus-storage split follows from the
demand and policy assumptions used in this model and does not generalise
to CCU deployment.

### Limitations affecting the financing-gap interpretation

The WTP benchmarks are imperfect proxies for target-scale demand. The
stated path uses a small purchaser survey and records intentions, which
may exceed realised commitments [@rodemeier_willingness_2023]. The
revealed path records signed deals but is sparse, concentrated among
early buyers, and includes first-of-a-kind premia. The 2030 gap is based
on early-market evidence. The 2040 and 2050 values are assumptions.

The benchmark simplifies demand-side heterogeneity. It applies one WTP
value across buyer sectors, CDR technologies, and target volumes, and it
averages BECCS and DACCS values equally. No study estimates a
price-dependent CDR demand curve for Europe, so the benchmark cannot be
built as one. While this matches the technology-neutral credit model, it
cannot capture method-specific preferences or declining willingness to
pay as volumes increase. The price-gap result therefore tests whether
benchmark prices cover the target-binding credit price. It does not test
whether enough buyers exist at that price to purchase the full target
volume.

The system-level gap is not an unfunded public subsidy requirement. Some
costs could be recovered through prices, tariffs, public infrastructure
finance, or contracts. The metric instead shows the additional
annualised system cost of reaching the target and its position in the
model. Country allocation assigns costs to the locations of investments
and operating changes in a least-cost European system. It does not
describe fiscal burden sharing or the final bill for governments,
consumers, or firms. This matters because the model treats Europe as
integrated, while policy authority is fragmented and major hosts,
especially the UK, sit outside EU burden-sharing and ETS governance.

### Future work

These limitations suggest several directions for future work. The most
immediate step is to model the EU ETS endogenously, including banking,
the Market Stability Reserve, and the integration of credited CDR into
the cap. Future research could also test alternative removal-target
trajectories and allocation rules, as the target level determines where
the supply curve is read and how strongly storage and system-cost
constraints apply.

The demand side can be improved by differentiating WTP across sectors
and between BECCS and DACCS, and by linking it to a volume-dependent
demand curve instead of a flat benchmark. The supply side can be refined
with project-level storage data, explicit transport-route constraints,
and country-specific financing assumptions. Where computing resources
permit, using perfect foresight, multiple weather years, finer temporal
resolution, and smaller credit-price steps would enhance timing, siting,
and threshold estimates.

# Conclusion {#sec:conclusion}

This study examined how CDR credit prices shape cost-optimal CDR
deployment in Europe and what financing gap remains. A standalone credit
price was added to the spatially resolved PyPSA-Eur model and swept
across low, medium, and high technology cost pathways. Each run records
the deployed volume, the technology mix, where capture and storage
locate, and the energy-system response. The resulting target-binding
prices were then compared with buyer willingness to pay, and the added
annualised system cost with the zero-credit baseline.

Target-scale CDR does not deploy without a dedicated CDR credit price.
In the zero-credit baseline, the ETS supports fossil & process CCS,
while atmospheric removal is absent or minimal. Introducing a credit
enables the CDR target to be reached at 100 to 150 €/t in both 2030,
2040 and 2050. Across all years and cost sensitivities, the
target-binding credit price remains within 100 /t--150 /t. Cost
sensitivities affect this price by no more than one 50 €/t increment. In
the scenarios tested, achieving the target depends more on the
credit-price threshold than on the technology-cost pathway. The
threshold is a modelled clearing price without risk premia or
certification costs, and is read as a floor.

Deployment shifts from BECCS to DACCS and concentrates in a few host
countries. The small 2030 CDR target is met entirely by BECCS within
available biomass. DACCS overtakes BECCS by 2040 in the medium
sensitivity and supplies 86 to 94% of removal by 2050 across all
sensitivities. Of the two DAC variants, only S-DAC is built, because the
higher electricity demand of L-DAC outweighs its lower capital cost at
European power prices. Denmark and the UK host 65% of capture in 2040,
and 70% of capture and 80% of storage in 2050. Storage access and
low-cost power set this geography, and the country-level DAC weather
factor shows no measurable effect on siting.

The main system effect is added capacity and annualised cost rather than
higher energy prices. DACCS at scale adds a large electricity, heat, and
biomass demand, with low-grade DAC heat the largest block. The model
meets this demand by building new wind and heat-pump capacity. European
average electricity prices rise by at most 5 €/MWh, and average heat
prices fall in 2050 as that capacity expands, while host-country prices
move more strongly in both directions. This stability assumes new
capacity can be built without friction.

Buyer willingness to pay covers the target-binding price in 2030 in all
but one corner case, low revealed WTP combined with high cost. Coverage
weakens as the target grows and the projected WTP paths diverge. Under
low WTP, a gap opens on the stated path from 2040 and reaches 4 to
54 €/t by 2050, while on the revealed path it appears earlier, widens to
26 €/t in 2040, and narrows to the high-cost case alone at 37 €/t by
2050. Under medium and high WTP, the price gap stays closed in every
year. The post-2030 WTP values are sensitivity assumptions, so the later
gaps carry the most uncertainty.

Meeting the target adds 46 bn €/yr of annualised system cost by 2050 in
the medium sensitivity, within 33 to 63 bn €/yr across cost
sensitivities. This cost concentrates where removal is hosted. The UK
carries about 22bn€/yr and Denmark about 12bn€/yr by 2050, driven by DAC
capital and transport and storage. Countries such as France pay close to
1 bn €/yr for renewable generation and biomass while hosting little
capture, and a few countries fall below the baseline as activity
relocates to the hubs. Credit revenue accrues to the capture operator,
so this wider cost is left to other instruments. About half of the
cost-optimal hosting sits in the UK, outside EU burden-sharing, which
makes the allocation question cross-border.

In the scenarios tested, a standalone carbon removal market delivers
cost-optimal deployment of CDR at Europe's target scale for credit
prices of 100 to 150 €/t. Early buyers cover that price under most
tested assumptions. Creating demand and allocating the added system cost
are separate problems. The credit resolves the first within the tested
range, while the second falls on a few host countries and several
non-hosting payers and needs a dedicated allocation mechanism.

Future work should make the EU ETS endogenous, test alternative target
trajectories, and differentiate willingness to pay and cost of capital
by buyer and country.

[]{#LastRealPage label="LastRealPage"}

# Appendix {#sec:appendix}

## Methodology {#subsec:app_supply_side}

This appendix documents the modelling detail behind the analysis,
organised to follow the three runs in the Results: the Baseline
reference system, Scenario 1 (deployment response to the credit price),
and Scenario 2 (the financing gap). The Baseline block describes the
base PyPSA-Eur model and the input data common to every run. The
Scenario 1 block covers the two model extensions together with the
techno-economic, resource, and price inputs that drive cost-optimal CDR
deployment across the credit-price sweep and the low, medium, and high
cost sensitivities. The Scenario 2 block builds the stated and revealed
willingness-to-pay benchmarks compared against the resulting credit
price. Figure [11](#fig:app_methodology_overview){reference-type="ref"
reference="fig:app_methodology_overview"} summarises this supply-side
and demand-side structure, which Scenario 2 integrates by setting the
target-binding credit price ($\mathrm{CP}^{*}$) and the system cost
against willingness to pay.
Table [7](#tab:core_assumptions){reference-type="ref"
reference="tab:core_assumptions"} lists the core assumptions and points
to the subsection documenting each.

<figure id="fig:app_methodology_overview" data-latex-placement="H">
<img src="./Figures/Appendix/Methodology_Overview.png" />
<figcaption>Overview of the supply-side and demand-side
methodology.</figcaption>
</figure>

::: {#tab:core_assumptions}
  **Assumption**            **Main values**                                              **Role in the analysis**                                                                       **App.**
  ------------------------- ------------------------------------------------------------ ---------------------------------------------------------------------------------------------- -------------------------------------------------------------------------------------------------------------------------------------
  Planning years            , 2040, 2050                                                 Years solved in the myopic PyPSA-Eur sequence.                                                 [7.1.1](#subsubsec:app_pypsa_config){reference-type="ref" reference="subsubsec:app_pypsa_config"}
  CDR credit grid           --500 EUR/tCO~2~, in 50 EUR/tCO~2~ steps                     Exogenous price sweep defining the CDR supply curve.                                           [7.1.3](#subsubsec:app_carbon_accounting){reference-type="ref" reference="subsubsec:app_carbon_accounting"}
  Cost sensitivities        Low, medium, high for BECCS, S-DAC and L-DAC                 Coherent cost, efficiency and energy-demand pathways.                                          [7.1.6](#subsubsec:app_tech_cost_scenarios){reference-type="ref" reference="subsubsec:app_tech_cost_scenarios"}
  Storage cost range        --20 EUR/tCO~2~, country-differentiated                      Spatial driver of storage and capture location.                                                [7.1.9](#subsubsec:app_sequestration_costs_potential){reference-type="ref" reference="subsubsec:app_sequestration_costs_potential"}
  CO~2~ storage potential   , 320, 600 MtCO~2~/yr                                        Annual infrastructure availability limit for geological sequestration.                         [7.1.9](#subsubsec:app_sequestration_costs_potential){reference-type="ref" reference="subsubsec:app_sequestration_costs_potential"}
  CDR policy target         , 160, 514 MtCO~2~/yr                                        Post-processed threshold used to identify $\mathrm{CP}^{*}$. Not an optimisation constraint.   [7.1.11](#subsubsec:app_policy_goals){reference-type="ref" reference="subsubsec:app_policy_goals"}
  EU ETS price              , 279, 463 EUR/tCO~2~                                        Exogenous carbon price for the wider system. Cancelled for eligible CDR revenue.               [7.1.12](#subsubsec:app_ets_prices){reference-type="ref" reference="subsubsec:app_ets_prices"}
  WTP benchmarks            Stated 197/250/302, revealed 136/381/612 EUR/tCO~2~ (2030)   Buyer-side benchmark for the price-gap analysis. 2040 and 2050 are projected sensitivities.    [7.1.13](#app:wtp-construction){reference-type="ref" reference="app:wtp-construction"}

  : Core assumptions used in the main analysis.
:::

**Baseline: Reference Model**\

------------------------------------------------------------------------

The Baseline is the cost-optimal European system at medium cost with no
CDR credit. The configuration and input data below define the base
PyPSA-Eur model used by all three runs.

### PyPSA-Eur model configuration {#subsubsec:app_pypsa_config}

<figure id="fig:pypsa_structure" data-latex-placement="H">
<img src="./Figures/Methodology/multisector_figure.png"
style="width:38.0%" />
<figcaption>Overview of sectors and links represented in
PyPSA-Eur</figcaption>
</figure>

#### Spatial and temporal resolution.

The model represents the European energy system at 96 nodes, clustered
from 35 countries by k-means on transmission topology and demand. This
resolution sits above the country level, so it resolves the
intra-country variation in renewable resources, biomass supply, and
geological storage cost that governs where CDR is sited, and it keeps
cross-border CO~2~ and electricity flows explicit. It stays coarse
enough to run the full sweep of credit prices, cost sensitivities, and
planning years. Capacity factors, ambient temperature, and heat demand
come from ERA5 reanalysis for 2019. The 8 760 hourly steps are reduced
to 168 representative segments with TSAM [@kotzur_impact_2018], each
weighted by the number of hours it represents. Transmission expansion is
capped at 1.5 times the cost-optimal volume.

#### Cost annualisation.

Investment costs are annualised over each technology's economic lifetime
using a uniform weighted average cost of capital of 7%, the PyPSA-Eur
default. The resulting annuity is added to the fixed and variable
operating costs to give the annualised system cost minimised by the
optimiser.

#### Sector coupling and energy carriers.

The model couples electricity, heating (district and building level),
transport (road, shipping, aviation), and industry. Energy carriers are
electricity, hydrogen, methane, liquid hydrocarbons, heat, and CO~2~.
Conversion technologies link these carriers endogenously across sectors:
electrolysis, methanation, fuel cells, Fischer-Tropsch synthesis, carbon
capture units, and direct air capture.
Figure [12](#fig:pypsa_structure){reference-type="ref"
reference="fig:pypsa_structure"} gives an overview of the represented
sectors and the conversion links between them.

#### CO~2~ bus hierarchy.

CO~2~ flows run through three carrier-level buses. The global
`co2 atmosphere` bus is the source for direct air capture and the sink
for biogenic re-emissions. Node-level `co2 stored` buses aggregate
captured CO~2~ from all sources (DAC, BECCS, point-source CCS) before
routing to pipelines or storage. Node-level `co2 sequestered` stores
represent permanent geological injection with country-specific capital
costs
(Appendix [7.1.9](#subsubsec:app_sequestration_costs_potential){reference-type="ref"
reference="subsubsec:app_sequestration_costs_potential"}).

#### Extendable technologies.

The following technologies may be invested in endogenously: onshore
wind, offshore wind (AC, DC, floating), solar PV (fixed-tilt and
horizontal single-axis tracking), solar rooftop, nuclear, open- and
combined-cycle gas turbines, biomass generators, batteries, hydrogen
stores, and all DACCS and BECCS capture links. Existing conventional
capacity (coal, lignite, hydro, geothermal) is fixed at reported 2020
installed levels and declines through decommissioning schedules. No
CO~2~ emissions cap is applied. CDR deployment is governed entirely by
the credit price signal and geological storage constraints.

#### Nuclear: scope and sensitivity note. {#subsubsec:app_nuclear_table}

Nuclear is an endogenously extendable technology, so the optimiser can
invest where it is cost-competitive, bounded by country-specific
capacity ceilings that reflect national policy and announced project
pipelines. The ceilings are aggregated capacity constraints on the
`Generator` component, tightened or relaxed by planning year to follow
scheduled phase-outs and new builds. Germany's fleet was decommissioned
in 2022 and is excluded. Belgium and Spain have legislated phase-outs
(2025 and 2027), so their ceilings are fixed at existing levels and the
remaining fleet is removed on schedule. All other nuclear-operating
countries keep existing capacity and may expand within the bounds in
Table [8](#tab:app_nuclear_limits){reference-type="ref"
reference="tab:app_nuclear_limits"}.

::: {#tab:app_nuclear_limits}
  --------------------- ---------------- ----------------------- ----------------------- -----------------------
  Country               Policy status                       2030                    2040                    2050
                                           \[GW$_{\mathrm{e}}$\]   \[GW$_{\mathrm{e}}$\]   \[GW$_{\mathrm{e}}$\]
  Belgium (BE)          Phase-out 2025                       4.1                     4.1                     4.1
  Bulgaria (BG)         Expansion                            2.1                     4.1                     6.1
  Switzerland (CH)      Existing fleet                       2.3                     2.3                     2.3
  Czech Republic (CZ)   Expansion                            5.4                     6.6                     8.6
  Spain (ES)            Phase-out 2027                       7.4                     7.4                     7.4
  Finland (FI)          Expansion                            5.8                     6.6                     7.4
  France (FR)           Expansion                           67.3                    75.0                    82.0
  United Kingdom (GB)   Expansion                           13.2                    16.4                    19.6
  Hungary (HU)          Expansion                            2.0                     4.4                     4.4
  Netherlands (NL)      Expansion                            0.5                     2.5                     4.5
  Romania (RO)          Expansion                            2.1                     3.5                     4.9
  Sweden (SE)           Expansion                            7.2                    10.0                    14.0
  Slovenia (SI)         Expansion                            0.7                     1.8                     1.8
  Slovakia (SK)         Existing fleet                       3.4                     3.4                     3.4
  Ukraine (UA)          Existing fleet                      13.8                    13.8                    13.8
  --------------------- ---------------- ----------------------- ----------------------- -----------------------

  : Country-specific upper bounds on nuclear installed capacity
  \[GW$_{\mathrm{e}}$\] for 2030, 2040, and 2050. Bounds reflect
  existing capacity plus announced project pipelines. Countries with a
  legislated phase-out (Belgium, Spain) are capped at existing installed
  levels. Germany is excluded from the model. No minimum capacity
  constraints are applied.
:::

#### Objective function.

Total annualised system cost is minimised:

$$\begin{equation}
\min \;
\sum_{i \in \mathcal{I}} \bigl(c_i^{\mathrm{inv}} + c_i^{\mathrm{fix}}\bigr)\,\bar{p}_i
\;+\;
\sum_{i \in \mathcal{I}} \sum_{t \in \mathcal{T}} w_t\, c_i^{\mathrm{var}}\, p_{i,t}
\;+\;
p^{\mathrm{CO_2}}\, E^{\mathrm{fossil}}
\;-\;
\sum_{o \in \mathcal{O}} p_o^{\mathrm{CDR}}\, C_o^{\mathrm{CDR}}
\label{eq:objective}
\end{equation}$$

where $\bar{p}_i$ is installed capacity of component $i$,
$c_i^{\mathrm{inv}}$ and $c_i^{\mathrm{fix}}$ are annualised capital and
fixed O&M costs (EUR MW$^{-1}$ yr$^{-1}$), $c_i^{\mathrm{var}}$ is the
variable cost (EUR MWh$^{-1}$), $p_{i,t}$ is dispatch at snapshot $t$,
$w_t$ is the snapshot weighting (hours represented by segment $t$),
$p^{\mathrm{CO_2}}$ is the exogenous EU ETS allowance price
(EUR tCO~2~$^{-1}$, see
Appendix [7.1.12](#subsubsec:app_ets_prices){reference-type="ref"
reference="subsubsec:app_ets_prices"}), $E^{\mathrm{fossil}}$ is total
annual fossil CO~2~ emissions, $p_o^{\mathrm{CDR}}$ is the CDR credit
price for origin $o$ (EUR tCO~2~$^{-1}$), $C_o^{\mathrm{CDR}}$ is the
annual credited CDR volume for origin $o$ (tCO~2~ yr$^{-1}$), and
$\mathcal{O} = \{\mathrm{dac},\, \mathrm{biogenic}\}$ is the set of
eligible CDR origins. The final term is implemented through the
accounting module described in
Appendix [7.1.3](#subsubsec:app_carbon_accounting){reference-type="ref"
reference="subsubsec:app_carbon_accounting"}.

### PyPSA-Eur default input data {#subsubsec:app_data_collection}

PyPSA-Eur provides a built-in, automated data pipeline that retrieves
and preprocesses core geographic, land-use, energy, and socio-economic
inputs on demand using a Snakemake-based workflow, with foundational
datasets supplied via a centralized Zenodo data bundle. Table
[9](#tab:auxiliary_datasets){reference-type="ref"
reference="tab:auxiliary_datasets"} gives an overview of auxiliary
datasets used in PyPSA-Eur.

::: {#tab:auxiliary_datasets}
  **Dataset type**    **Content**                         **Usage**
  ------------------- ----------------------------------- ----------------------------------------
  Geographic          NUTS3 shapes, EEZ boundaries        Administrative regions, maritime zones
  Land use            CORINE land cover, Natura 2000      Renewable energy site eligibility
  Energy statistics   Historical hydropower generation    Hydropower modelling validation
  Demographics        GDP and population at NUTS3 level   Spatial demand distribution

  : Overview of auxiliary datasets used in PyPSA-Eur
:::

**Scenario 1: Deployment Response**\

------------------------------------------------------------------------

Scenario 1 introduces the standalone CDR credit and sweeps its price to
trace the supply curve. The two model extensions enter here: Module 1,
the standalone credit market
(Appendix [7.1.3](#subsubsec:app_carbon_accounting){reference-type="ref"
reference="subsubsec:app_carbon_accounting"}), and Module 2, the
spatially resolved CDR representation, which spans the DACCS variants
(Appendix [7.1.5](#subsubsec:app_daccs_variants){reference-type="ref"
reference="subsubsec:app_daccs_variants"}), the weather-dependent DAC
energy demand (Appendix [7.1.7](#app:dac-weather){reference-type="ref"
reference="app:dac-weather"}), and the country-differentiated storage
costs
(Appendix [7.1.9](#subsubsec:app_sequestration_costs_potential){reference-type="ref"
reference="subsubsec:app_sequestration_costs_potential"}). The remaining
blocks set the cost sensitivity, resource limits, CDR target, and carbon
price.

### Standalone CDR Credit Market (Module 1) {#subsubsec:app_carbon_accounting}

#### Problem statement.

Standard PyPSA-Eur does not distinguish CO~2~ by origin. Without
modification, any CO~2~ delivered to a `co2 sequestered` store generates
the same cost reduction, meaning fossil carbon capture would earn CDR
credit revenue indistinguishably from atmospheric or biogenic sources.
The module addresses this by restricting credits to eligible origins and
enforcing annual market demand limits within the solver.

#### Origin classification.

Each capture link is assigned an origin from its PyPSA-Eur carrier name
using a deterministic string classifier. Links whose carrier contains
`dac` are classified as `dac`. Links containing any of the tokens
`biomass`, `biogas`, `biosng`, `btl`, `fuelwood`, `msw`, or `bio` are
classified as `biogenic`. All remaining capture links are classified as
`fossil`. CO~2~ infrastructure links (pipelines, sequestration transfer
links) are excluded. The credit-eligible set is controlled by the
`cdr_credit_scope` configuration parameter, set to `["dac", "biogenic"]`
in all supply curve scenarios, explicitly excluding fossil carbon
capture. These two eligible origins are the set
$\mathcal{E} = \{\text{BECCS},\ \text{DACCS}\}$ of
Eq. [\[eq:eligible\]](#eq:eligible){reference-type="ref"
reference="eq:eligible"}, where DACCS corresponds to the `dac` origin
and BECCS is used as shorthand for the full `biogenic` capture bucket,
which also covers biogas upgrading, BioSNG, biomass-to-liquid, and
municipal-waste capture.

#### Decision variables.

For each eligible origin $o \in \mathcal{E}$, two scalar decision
variables are added at solve time:

- $s^{\mathrm{seq}}_{o}$: annual CO~2~ attributed to geological
  sequestration for origin $o$ (tCO~2~ yr$^{-1}$)

- $s^{\mathrm{cred}}_{o}$: annual credited CDR volume for origin $o$,
  bounded above by the market demand cap (tCO~2~ yr$^{-1}$)

#### Constraints.

Three constraint sets enforce physical consistency:

$$\begin{align}
  \sum_{o \in \mathcal{E}} s^{\mathrm{seq}}_{o} &\leq S^{\mathrm{seq}}
  \label{eq:cdr_total_seq} \\
  s^{\mathrm{seq}}_{o} &\leq K_{o} \quad \forall\, o \in \mathcal{E}
  \label{eq:cdr_per_origin} \\
  s^{\mathrm{cred}}_{o} &\leq \min\!\bigl(s^{\mathrm{seq}}_{o},\; K_{o}\bigr)
    \quad \forall\, o \in \mathcal{E}
  \label{eq:cdr_credited}
\end{align}$$

where
$S^{\mathrm{seq}} = \sum_t w_t \sum_{l \in \mathcal{L}^{\mathrm{seq}}} p_{l,t}$
is total physical annual sequestration across all `co2 sequestered`
links $\mathcal{L}^{\mathrm{seq}}$, and
$K_o = \sum_t w_t \sum_{i:\,\mathrm{origin}(i)=o} \eta_i^{\mathrm{CO_2}}\, p_{i,t}$
is annual CO~2~ capture by origin $o$, with $\eta_i^{\mathrm{CO_2}}$ the
link's CO~2~ output coefficient and $w_t$ the snapshot weighting.
Eqs. [\[eq:cdr_per_origin\]](#eq:cdr_per_origin){reference-type="ref"
reference="eq:cdr_per_origin"} and
[\[eq:cdr_credited\]](#eq:cdr_credited){reference-type="ref"
reference="eq:cdr_credited"} are the solve-time form of the bounds
$s^{\mathrm{cred}}_{o} \le K_{o}$ and
$\sum_{o} s^{\mathrm{cred}}_{o} \le S^{\mathrm{seq}}$ stated in
Eq. [\[eq:eligible\]](#eq:eligible){reference-type="ref"
reference="eq:eligible"}. An additional constraint caps total credited
volume against the annual market demand limit $D^{\mathrm{cap}}$:
$\sum_{o \in \mathcal{E}} s^{\mathrm{cred}}_{o} \leq D^{\mathrm{cap}}$.
In the supply curve runs $D^{\mathrm{cap}}$ is set non-binding, so the
sweep traces an unconstrained price response; the constraint is retained
in the formulation for scenarios with an explicit demand trajectory.

#### Objective adjustment.

CDR credit revenue is subtracted from the system cost objective at solve
time:

$$\begin{equation}
  \Delta z \;=\; -\,\mathrm{CP} \sum_{o \in \mathcal{E}} s^{\mathrm{cred}}_{o}
\end{equation}$$

This is the credit term $C^{\mathrm{credit}}$ of
Eq. [\[eq:credit\]](#eq:credit){reference-type="ref"
reference="eq:credit"} and corresponds to the final term in
Equation [\[eq:objective\]](#eq:objective){reference-type="ref"
reference="eq:objective"}. The implementation also supports
origin-differentiated credit prices, but all runs in this study apply a
single uniform $\mathrm{CP}$.

#### ETS interaction.

With `cdr_credit_standalone: true`, DAC and BECCS operators are treated
as outside the EU ETS and do not hold allowances. In the default
PyPSA-Eur formulation, atmospheric CO~2~ withdrawal would implicitly
generate an ETS benefit through the shadow value of the emissions
constraint. To cancel this double-counting, a marginal cost surcharge of
$p^{\mathrm{ETS}}_{t} \cdot \eta_i^{\mathrm{atm}}$ is applied to each
eligible capture link at network construction time, where
$\eta_i^{\mathrm{atm}}$ is the atmospheric withdrawal coefficient per
unit of dispatch. Summed over nodes and periods, these surcharges equal
the cancellation cost $C^{\mathrm{cancel}}$ of
Eq. [\[eq:cancel\]](#eq:cancel){reference-type="ref"
reference="eq:cancel"}, with $a_{n,t}$ the aggregate of
$\eta_i^{\mathrm{atm}}\, p_{i,t}$ over the eligible links at node $n$.
CO~2~ routed to carbon capture and utilisation retains ETS neutrality:
the implicit credit at capture is offset by the ETS liability at
re-emission.

### Credited-CDR attribution (waterfall reconciliation) {#subsubsec:app_waterfall}

The solve-time module of
Appendix [7.1.3](#subsubsec:app_carbon_accounting){reference-type="ref"
reference="subsubsec:app_carbon_accounting"} lets the credit act on
eligible removals through the objective, so BECCS and DAC deployment
respond to the credit price. What that module leaves open is how the
credited total divides between the two eligible origins. In PyPSA-Eur
all captured , whether fossil, biogenic, or DAC, is commingled in a
single `co2 stored` pool before part is sequestered and part is drawn
back out for synthetic fuels. Once pooled there is no physical fact
about which captured tonne ends up underground, so under a single credit
price the solver is indifferent to how the credited volume is split
between DAC and biogenic: the per-origin split is degenerate. It is
fixed after solving by a deterministic reconciliation, applied to the
solved flows without re-optimisation, that reproduces the incentive
already priced into the objective.

Eligible (biogenic and DAC) removals take first claim on the physical
geological sequestration of the year; fossil is stored only with the
residual capacity and earns no credit. Let $A_{\mathrm{dac}}$,
$A_{\mathrm{bio}}$, and $A_{\mathrm{fos}}$ be the annual captured by
each origin, $A_{\mathrm{elig}} = A_{\mathrm{dac}} +
A_{\mathrm{bio}}$ the eligible capture, and $Q^{\mathrm{seq}}$ the total
physical sequestration. The credited eligible removal is
$$\begin{equation}
    s^{\mathrm{cred}} = \min\!\left(A_{\mathrm{elig}},\, Q^{\mathrm{seq}}\right),
    \qquad
    C_o^{\mathrm{CDR}} = s^{\mathrm{cred}}\,\frac{A_{o}}{A_{\mathrm{elig}}}
    \quad \text{for } o \in \{\mathrm{dac},\, \mathrm{bio}\},
    \label{eq:waterfall}
\end{equation}$$ with the fossil residual
$\min\!\left(A_{\mathrm{fos}},\, Q^{\mathrm{seq}} -
s^{\mathrm{cred}}\right)$ stored but uncredited. When eligible capture
exceeds available storage, the scarce storage is shared within the
eligible bucket in proportion to captured volume. By construction the
rule satisfies $$\begin{equation}
    s^{\mathrm{cred}} \le A_{\mathrm{elig}}
    \qquad\text{and}\qquad
    s^{\mathrm{cred}} \le Q^{\mathrm{seq}},
    \label{eq:waterfall_bounds}
\end{equation}$$ so credited CDR can never exceed either eligible
capture or what is physically stored. These are the same bounds the
optimiser enforces at solve time
(Eqs. [\[eq:cdr_total_seq\]](#eq:cdr_total_seq){reference-type="ref"
reference="eq:cdr_total_seq"} to
[\[eq:cdr_credited\]](#eq:cdr_credited){reference-type="ref"
reference="eq:cdr_credited"}), so the reconciliation reads off a
quantity the solver already valued rather than imposing an external
assumption: the credit lets eligible removals out-bid fossil for scarce
storage, so $s^{\mathrm{cred}}$ equals the eligible volume priced into
deployment. The credited total
$s^{\mathrm{cred}} = \sum_{o} C_o^{\mathrm{CDR}}$ is the quantity
compared against the yearly target to read off the target-binding credit
price $\mathrm{CP}^{*}$
(Eq. [\[eq:cp_star\]](#eq:cp_star){reference-type="ref"
reference="eq:cp_star"}).

### DACCS technology variants {#subsubsec:app_daccs_variants}

This work models two DACCS technology variants reflecting distinct
process routes for atmospheric CO$_2$ capture. The variants differ in
sorbent type, energy carrier requirements, and the degree to which they
produce recoverable heat, leading to different cost structures and
synergies with the broader energy system.

**Solid-sorbent DAC (S-DAC)** employs temperature--vacuum swing
adsorption (TVSA) on solid amine or zeolite sorbents. Desorption is
driven by low-grade heat (approximately 80--120 °C), making S-DAC
well-suited for integration with district heating networks or industrial
waste heat. Electricity consumption is modest. The process outputs
recoverable heat that can be fed back into local heat buses.

**Liquid-sorbent DAC (L-DAC)** circulates an aqueous KOH solution to
absorb atmospheric CO$_2$, which is subsequently liberated through
high-temperature calcination. In the all-electric configuration modelled
here, the calciner is powered electrically, resulting in substantially
higher electricity consumption than S-DAC but no requirement for a
separate low-temperature heat supply. A modest amount of waste heat is
still recovered.

Beyond the capture step, both variants require compression to bring
captured CO$_2$ to pipeline pressure, adding
0.11--0.13 MWh$_{\mathrm{el}}$/tCO$_2$ of electricity and
0.16--0.19 MWh/tCO$_2$ of recoverable compression heat across variants.

In PyPSA-Eur each variant is added as a multi-bus `Link` via the
`add_dac` routine, connecting the regional electricity bus (`bus0`), the
local heat bus (`bus1`), the atmospheric CO$_2$ bus (`bus2`, source),
and the regional CO$_2$ storage bus (`bus3`, sink). Link efficiencies
are derived per node from the base energy inputs together with
country-level weather scaling factors (applied to S-DAC and L-DAC, see
Appendix [7.1.7](#app:dac-weather){reference-type="ref"
reference="app:dac-weather"}). Capital costs are normalised to
EUR MW$^{-1}$ of electricity input so that the PyPSA-Eur investment
framework handles capacity correctly. Recovered S-DAC heat (about 50 °C)
is credited to the local heat bus at 0.52 of its full district-heating
value, an exergy-quality adjustment that reflects its low temperature.
The optimisation selects endogenously among variants based on their
levelised cost and regional energy conditions. The techno-economic
parameters for both variants across cost sensitivities and planning
years are reported together with the other capture technologies in
Appendix [7.1.6](#subsubsec:app_tech_cost_scenarios){reference-type="ref"
reference="subsubsec:app_tech_cost_scenarios"}.

### CDR technology cost sensitivities {#subsubsec:app_tech_cost_scenarios}

This appendix documents the construction of the low, medium, and high
technology cost sensitivities referenced in
Section [3.2](#subsec:scenarios){reference-type="ref"
reference="subsec:scenarios"} and applied to the three carbon capture
technologies represented in the model: biomass CHP capture (BECCS),
solid-sorbent DAC (S-DAC), and liquid-sorbent DAC (L-DAC).

#### Data sources.

The medium sensitivity reflects the medium estimate of the Danish Energy
Agency (DEA) technology catalogue for carbon capture, transport, and
storage [@danish_energy_agency_technology_2024], which is the same data
source used by default in PyPSA-Eur. The low and high sensitivtiies for
2050 are also drawn directly from the DEA catalogue, which reports a low
and high cost range for the 2050 horizon for each technology. For 2030
and 2040, the DEA catalogue does not report low and high sensitivtiies
for every parameter. These are constructed using the percentage-spread
method described below.

#### Construction of low and high sensitivtiies.

For each technology and parameter, the low and high values for 2050 are
taken directly from the DEA technology catalogue. For 2030 and 2040, the
same relative spread observed in 2050 is applied to the medium scenario
values. For a parameter $x$ in year $y$, the sensitivtiies are
constructed as

$$\begin{equation}
x_{\text{low}}^{y} = x_{\text{med}}^{y} \cdot
   \left(1 - \frac{x_{\text{med}}^{2050} - x_{\text{low}}^{2050}}{x_{\text{med}}^{2050}}\right),
\qquad
x_{\text{high}}^{y} = x_{\text{med}}^{y} \cdot
   \left(1 + \frac{x_{\text{high}}^{2050} - x_{\text{med}}^{2050}}{x_{\text{med}}^{2050}}\right)
\label{eq:tech_cost_envelope}
\end{equation}$$

for $y \in \{2030, 2040\}$. This preserves the proportional uncertainty
range reported by the DEA while reflecting the generally higher costs
and energy demands of earlier deployment years.

The technologies differ in investment cost, fixed O&M (FOM), electricity
and heat demand, heat recovery, and lifetime. Compression energy for
CO$_2$ transport and storage is modelled separately.

#### Compression parameters.

The current DEA catalogue revision does not report scenario- or
year-specific values for compression electricity demand and recoverable
compression heat. These parameters are therefore taken from the previous
DEA catalogue edition and held constant across years.

For biomass CHP capture, the previous catalogue reported a low, medium,
and high range, which is retained here. For the DACCS variants, only a
single medium value was available and is therefore applied uniformly
across all scenarios. As compression energy represents only a small
share of total DACCS energy demand, the resulting sensitivity is
limited.

#### Parameter tables.

Tables [10](#tab:app_tech_cost_beccs){reference-type="ref"
reference="tab:app_tech_cost_beccs"} to
[12](#tab:app_tech_cost_ldac){reference-type="ref"
reference="tab:app_tech_cost_ldac"} summarise the resulting
techno-economic assumptions by technology, year, and scenario.
Investment costs are reported in EUR/(tCO$_2$/h), FOM as annual
percentage of investment cost, energy inputs and outputs in MWh/tCO$_2$,
and lifetime in years.

::: {#tab:app_tech_cost_beccs}
+--------------------+-----------------------------------+-----------------------------------+-----------------------------------+
|                    | low                               | medium                            | high                              |
+:===================+:=========:+:=========:+:=========:+:=========:+:=========:+:=========:+:=========:+:=========:+:=========:+
| 2-4(lr)5-7(lr)8-10 | 2030      | 2040      | 2050      | 2030      | 2040      | 2050      | 2030      | 2040      | 2050      |
| Parameter          |           |           |           |           |           |           |           |           |           |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| FOM \[%/yr\]       | 2.87      | 2.86      | 2.86      | 3.02      | 3.00      | 3.00      | 3.04      | 3.03      | 3.03      |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Capture rate \[-\] | 0.85      | 0.90      | 0.90      | 0.90      | 0.95      | 0.95      | 0.90      | 0.95      | 0.95      |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Comp. elec. in     | 0.079     | 0.070     | 0.070     | 0.085     | 0.075     | 0.075     | 0.102     | 0.090     | 0.090     |
| \[MWh/t\]^a^       |           |           |           |           |           |           |           |           |           |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Comp. heat out     | 0.129     | 0.120     | 0.120     | 0.140     | 0.130     | 0.130     | 0.162     | 0.150     | 0.150     |
| \[MWh/t\]^a^       |           |           |           |           |           |           |           |           |           |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Elec. input        | 0.018     | 0.017     | 0.017     | 0.025     | 0.023     | 0.023     | 0.027     | 0.025     | 0.025     |
| \[MWh/t\]          |           |           |           |           |           |           |           |           |           |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Heat input         | 0.60      | 0.55      | 0.55      | 0.72      | 0.66      | 0.66      | 0.79      | 0.72      | 0.72      |
| \[MWh/t\]          |           |           |           |           |           |           |           |           |           |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Heat output        | 0.60      | 0.55      | 0.55      | 0.72      | 0.66      | 0.66      | 0.79      | 0.72      | 0.72      |
| \[MWh/t\]          |           |           |           |           |           |           |           |           |           |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Investment         | 2.34      | 2.10      | 2.10      | 6.30      | 5.66      | 5.66      | 8.46      | 7.60      | 7.60      |
| \[M€/(t/h)\]       |           |           |           |           |           |           |           |           |           |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| Lifetime \[yr\]    | 25        | 25        | 25        | 25        | 25        | 25        | 25        | 25        | 25        |
+--------------------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
| ^a^From previous DEA catalogue edition.                                                                                        |
+--------------------------------------------------------------------------------------------------------------------------------+

: Biomass CHP capture (BECCS): techno-economic parameters by scenario
and year. Compression-electricity-input and compression-heat-output
values are taken from the previous DEA catalogue edition.
:::

::: {#tab:app_tech_cost_sdac}
+--------------------+-----------------------+-----------------------+-----------------------+
|                    | low                   | medium                | high                  |
+:===================+:=====:+:=====:+:=====:+:=====:+:=====:+:=====:+:=====:+:=====:+:=====:+
| 2-4(lr)5-7(lr)8-10 | 2030  | 2040  | 2050  | 2030  | 2040  | 2050  | 2030  | 2040  | 2050  |
| Parameter          |       |       |       |       |       |       |       |       |       |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| FOM \[%/yr\]       | 3.61  | 3.68  | 3.87  | 3.92  | 4.00  | 4.21  | 3.78  | 3.85  | 4.05  |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| Comp. elec. in     | 0.11  | 0.11  | 0.11  | 0.11  | 0.11  | 0.11  | 0.11  | 0.11  | 0.11  |
| \[MWh/t\]^a^       |       |       |       |       |       |       |       |       |       |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| Comp. heat out     | 0.16  | 0.16  | 0.16  | 0.16  | 0.16  | 0.16  | 0.16  | 0.16  | 0.16  |
| \[MWh/t\]^a^       |       |       |       |       |       |       |       |       |       |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| Elec. input        | 0.166 | 0.144 | 0.130 | 0.230 | 0.200 | 0.180 | 0.345 | 0.300 | 0.270 |
| \[MWh/t\]          |       |       |       |       |       |       |       |       |       |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| Heat input         | 1.05  | 0.89  | 0.77  | 1.50  | 1.29  | 1.10  | 2.25  | 1.90  | 1.65  |
| \[MWh/t\]          |       |       |       |       |       |       |       |       |       |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| Heat output        | 0.93  | 0.70  | 0.70  | 1.00  | 0.75  | 0.75  | 1.33  | 1.00  | 1.00  |
| \[MWh/t\]          |       |       |       |       |       |       |       |       |       |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| Investment         | 2.59  | 1.81  | 1.55  | 3.57  | 2.50  | 2.14  | 5.36  | 3.75  | 3.21  |
| \[M€/(t/h)\]       |       |       |       |       |       |       |       |       |       |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| Lifetime \[yr\]    | 25    | 30    | 30    | 25    | 30    | 30    | 25    | 30    | 30    |
+--------------------+-------+-------+-------+-------+-------+-------+-------+-------+-------+
| ^a^Single medium value applied across all scenarios.                                       |
+--------------------------------------------------------------------------------------------+

: Solid-sorbent DAC (S-DAC): techno-economic parameters by scenario and
year. Compression parameters are held constant across scenarios.
:::

::: {#tab:app_tech_cost_ldac}
+--------------------+--------------------+--------------------+--------------------+
|                    | low                | medium             | high               |
+:===================+:====:+:====:+:====:+:====:+:====:+:====:+:====:+:====:+:====:+
| 2-4(lr)5-7(lr)8-10 | 2030 | 2040 | 2050 | 2030 | 2040 | 2050 | 2030 | 2040 | 2050 |
| Parameter          |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| FOM \[%/yr\]       | 3.79 | 4.03 | 3.73 | 3.65 | 3.88 | 3.59 | 3.81 | 4.05 | 3.75 |
+--------------------+------+------+------+------+------+------+------+------+------+
| Comp. elec. in     | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 | 0.13 |
| \[MWh/t\]^a^       |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Comp. heat out     | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 | 0.19 |
| \[MWh/t\]^a^       |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Elec. input        | 0.97 | 0.97 | 0.92 | 1.39 | 1.39 | 1.32 | 2.07 | 2.07 | 1.97 |
| \[MWh/t\]          |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Heat input         | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| \[MWh/t\]          |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Heat output        | 0.10 | 0.10 | 0.10 | 0.30 | 0.30 | 0.30 | 0.50 | 0.50 | 0.50 |
| \[MWh/t\]          |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Investment         | 2.26 | 1.59 | 1.34 | 3.29 | 2.32 | 1.95 | 4.94 | 3.49 | 2.93 |
| \[M€/(t/h)\]       |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Lifetime \[yr\]    | 25   | 25   | 25   | 30   | 30   | 30   | 35   | 35   | 35   |
+--------------------+------+------+------+------+------+------+------+------+------+
| ^a^Single medium value applied across all scenarios.                              |
+-----------------------------------------------------------------------------------+

: Liquid-sorbent DAC (L-DAC): techno-economic parameters by scenario and
year. Compression parameters are held constant across scenarios.
:::

#### Data quality notes and acknowledged inconsistencies.

A small number of inconsistencies in the DEA source data are retained
for transparency:

- **S-DAC FOM:** The DEA catalogue reports a 2050 high value below the
  medium value. Applying
  Equation [\[eq:tech_cost_envelope\]](#eq:tech_cost_envelope){reference-type="ref"
  reference="eq:tech_cost_envelope"} therefore results in slightly lower
  high values for 2030 and 2040. The original DEA values are retained
  unchanged, as the effect on total system cost is negligible.

- **L-DAC FOM:** Both the low and high 2050 values reported by the DEA
  are above the medium value, meaning the sensitivtiies do not fully
  bracket the medium sensitivity. The values are preserved without
  manual adjustment.

- **Biomass CHP capture rate:** The DEA catalogue reports the same
  medium and high capture rate in 2050 (95 %), resulting in no spread
  for the high scenario. This likely reflects an assumed technical upper
  bound for post-combustion capture efficiency.

These inconsistencies originate from the DEA catalogue itself rather
than from the envelope construction method applied here.

### Weather-dependent DAC energy demand {#app:dac-weather}

The energy demand of Direct Air Capture (DAC) depends on ambient
temperature and relative humidity (RH), which affect adsorption and
desorption processes. Following the methodology of
 [@okosun_global_2026], this relationship is represented by expressing
DAC energy requirements as a function of $(T,\mathrm{RH})$, with RH
derived from temperature and dew point data.

Temperature and dew point data from ERA5 (2 m height) are used to
compute RH using the Sonntag formula for saturation vapour pressure, as
in [@okosun_global_2026]:

$$\begin{equation}
  \ln e_\mathrm{w}(T) =
    \frac{-6096.9385}{T}
    - 2.711193 \times 10^{-2} T
    + 1.673952 \times 10^{-5} T^{2}
    + 2.433502 \ln T
    + 21.2409642,
\end{equation}$$

$$\begin{equation}
  \mathrm{RH} = 100 \times
    \exp\!\bigl[\ln e_\mathrm{w}(T_d) - \ln e_\mathrm{w}(T)\bigr],
\end{equation}$$

with RH clipped to $[0,100]\%$.

Country-level scaling factors are constructed by evaluating the energy
function $E(T,\mathrm{RH})$ at annual average conditions and normalising
to reference values ($T_{\mathrm{ref}} = 15^\circ\mathrm{C}$,
$\mathrm{RH}_{\mathrm{ref}} = 60\%$), consistent with
[@okosun_global_2026]:

$$\begin{equation}
  f^c = \frac{E(T_c, \mathrm{RH}_c)}{E(T_{\mathrm{ref}}, \mathrm{RH}_{\mathrm{ref}})}.
\end{equation}$$

Scaling factors are limited to $[0.7, 1.3]$ to avoid extrapolation
beyond the range covered in the source study.

Node-level energy demand is then given by:

$$\begin{equation}
  E^n = E^0 \cdot f^{c(n)},
\end{equation}$$

where each node $n$ inherits the factor of its corresponding country.

In PyPSA-Eur, DAC is implemented as a link converting electricity into
stored CO$_2$. The weather-dependent energy demand is incorporated by
adjusting link efficiencies according to $E^n$, with capital costs
scaled consistently.

### Biomass supply and non-biogenic capture build-rate limits {#subsubsec:app_biomass_constraints}

Solid biomass and biogas potentials are drawn from the JRC ENSPRESO
database [@ENSPRESO] at NUTS2 resolution and disaggregated to model
nodes. Each resource stream is a PyPSA-Eur generator whose annual
throughput is bounded by an `e_sum_max` constraint set to the regional
potential, taken here at the full ENSPRESO central estimate. The share
classified as unsustainable (e.g. primary forest harvesting) is phased
out linearly between 2025 and 2040 and fixed by a binding equality
constraint, so the model cannot avoid retiring it.

BECCS itself is not capped. Its deployment is bounded indirectly by this
biomass supply ceiling: once biomass is fully allocated across heat,
CHP, and industry, no further BECCS throughput is possible regardless of
the credit price. Build-rate limits (`limit_max_growth`) apply only to
the fossil and industrial capture carriers: process emissions CC, SMR
CC, gas CHP CC, and gas for industry CC.

### Spatially differentiated CO~2~ sequestration costs and potential {#subsubsec:app_sequestration_costs_potential}

#### Implementation in PyPSA-Eur.

Geological CO~2~ sequestration is represented through *co2 sequestered*
stores at each model node, connected to upstream *co2 stored* buses via
unit-efficiency links. Country-specific storage costs are implemented as
annualised capital costs per tonne of CO~2~ stored. These costs are
interpreted as levelised storage costs assuming a 50-year site lifetime
and a 7 % discount rate. Storage capacity at each node is constrained by
geological availability, and total system-wide sequestration is limited
through the global storage potential in
Table [13](#tab:co2-sequestration-potential){reference-type="ref"
reference="tab:co2-sequestration-potential"} and the crediting
constraints described in
Appendix [7.1.3](#subsubsec:app_carbon_accounting){reference-type="ref"
reference="subsubsec:app_carbon_accounting"}.

#### Cost typology and source basis.

CO~2~ storage costs are assigned based on geological formation type and
regional location, following the cost typology established by the Zero
Emissions Platform [@zero_emissions_platform_costs_2011] and
corroborated by the European Commission Joint Research Centre
[@martinez_castilla_clean_2025]. Both sources report storage costs in
the range of approximately 1--20 € / tCO$_2$, varying primarily with
three factors: (i) reservoir type (depleted oil and gas fields, DOGF,
versus saline aquifers, SA), (ii) onshore versus offshore location, and
(iii) reservoir characteristics such as injectivity and the availability
of legacy well infrastructure for reuse. Offshore saline aquifers
without legacy infrastructure represent the upper end of the range,
while offshore depleted fields with existing well infrastructure
represent the lower end. A JRC systematic review independently confirms
the same offshore range, 2--20 € / tCO$_2$ depending on formation type
and infrastructure availability [@JRC_CETO_2024].

Since CO~2~ transport is modelled separately as a pipeline network
(Appendix [7.1.10](#subsubsec:app_co2_transport){reference-type="ref"
reference="subsubsec:app_co2_transport"}), the values in
Table [\[tab:co2-seq-costs\]](#tab:co2-seq-costs){reference-type="ref"
reference="tab:co2-seq-costs"} represent pure geological storage costs
only, excluding compression and transport to the storage site.

#### Country-level mapping.

Country-level storage costs are assigned using a simplified mapping
approach. Each country is linked to a representative offshore region
(e.g. North Sea, Baltic Sea, Mediterranean, Adriatic, Black Sea,
Atlantic) and a formation type (DOGF or SA), as indicated in the
*Formation Type* column of
Table [\[tab:co2-seq-costs\]](#tab:co2-seq-costs){reference-type="ref"
reference="tab:co2-seq-costs"}. Countries with access to mature storage
infrastructure in the North Sea, characterised by abundant depleted oil
and gas fields and established injection experience, are assigned costs
at the lower end of the reported range. Countries relying on less
developed regions or saline aquifer storage are assigned correspondingly
higher values within the range. These assignments should be interpreted
as stylised cost assumptions rather than precise country-specific
estimates, reflecting the limited availability of site-level data for
most European countries.

For countries without direct access to offshore storage, a fallback
value of 20 € / tCO$_2$ is applied, representing the upper bound of the
reported range and implicitly reflecting the geological and logistical
constraints associated with onshore or landlocked storage contexts.

::: tabular
\@ll S\[table-format=2.0\] @ ll S\[table-format=2.0\] @ Country &
Formation & Cost \[€/t\] & Country & Formation & Cost \[€/t\]\
Norway & North Sea, DOGF (legacy) & 5 & Poland & Baltic, SA & 11\
Denmark & North Sea, DOGF/SA hub & 6 & Portugal & Atlantic, SA & 11\
Netherlands & North Sea, DOGF hub & 6 & Croatia & Adriatic, SA & 11\
United Kingdom & North Sea, DOGF (legacy) & 6 & Albania & Adriatic, SA &
12\
Belgium & North Sea, DOGF & 7 & Finland & Baltic, SA & 12\
Germany & North Sea, DOGF & 8 & Romania & Black Sea, SA & 12\
Ireland & Atlantic, SA & 8 & Montenegro & Adriatic, SA & 12\
France & Atlantic, SA & 9 & Ukraine & Black Sea, SA & 12\
Spain & Atlantic/Med., SA & 10 & Bulgaria & Black Sea, SA & 13\
Italy & Mediterranean, SA & 10 & Estonia & Baltic, SA & 13\
Sweden & Baltic/North Sea & 10 & Lithuania & Baltic, SA & 13\
Greece & Mediterranean, SA & 10 & Latvia & Baltic, SA & 13\
:::

*All other countries (no direct offshore access): fallback 20 €/tCO~2~.*

#### CO~2~ sequestration potential.

The availability of CO~2~ sequestration is represented through an
exogenously defined system-wide sequestration potential, which limits
the maximum amount of CO~2~ that can be stored annually. The assumed
trajectory reflects a gradual scale-up of CO~2~ transport and storage
infrastructure in Europe, consistent with policy targets outlined in the
European Commission's CCUS strategy. The annual sequestration potentials
applied in the model are listed in
Table [13](#tab:co2-sequestration-potential){reference-type="ref"
reference="tab:co2-sequestration-potential"}. These values act as upper
bounds on total annual CO~2~ sequestration in the model and represent
infrastructure availability constraints rather than geological storage
limits.

::: {#tab:co2-sequestration-potential}
  **Year**     **Potential (MtCO~2~/yr)**
  ---------- ----------------------------
  2025                                 10
  2030                                 70
  2035                                170
  2040                                320
  2045                                450
  2050                                600

  : Assumed CO~2~ sequestration potential over time.
:::

### CO~2~ transport network {#subsubsec:app_co2_transport}

Captured CO~2~ need not be stored where it is captured. With
`co2_spatial` and `co2_network` enabled, each node carries its own
`co2 stored` bus, and these buses are joined by an extendable,
bidirectional pipeline network that routes CO~2~ from capture sites to
storage-rich nodes. This routing lets DAC and BECCS locate apart from
the geological stores they feed, and it sets the
transport-versus-storage trade-off behind the spatial results in
Appendix [7.2.2](#app:scenario_tables){reference-type="ref"
reference="app:scenario_tables"}.

Each link connects two neighbouring `co2 stored` buses and is sized
endogenously (`p_nom_extendable`, with `p_min_pu` $=-1$ for
bidirectional flow). Its capital cost scales with length and with the
underwater fraction of the route, blending onshore and submarine
pipeline costs:

$$\begin{equation}
c^{\mathrm{pipe}}_{\ell} = f_y\, L_{\ell} \left[(1-u_{\ell})\,c^{\mathrm{on}} + u_{\ell}\,c^{\mathrm{sub}}\right],
\label{eq:co2_pipeline_cost}
\end{equation}$$

where $L_{\ell}$ is link length, $u_{\ell}$ the underwater fraction,
$c^{\mathrm{on}}$ and $c^{\mathrm{sub}}$ the onshore and submarine unit
investment costs (Table [14](#tab:app_co2_pipeline){reference-type="ref"
reference="tab:app_co2_pipeline"}), and $f_y$ a year-dependent network
cost factor. The factor falls from 5.0 in 2030 to 3.0 in 2040 and 1.5 in
2050, reflecting a declining build-out premium as the European CO~2~
transport network matures.

::: {#tab:app_co2_pipeline}
+----------------------------+--------------------+--------------------+
| Parameter                  | Onshore            | Submarine          |
+:===========================+:==================:+:==================:+
| Investment                 | 2116               | 4233               |
| \[EUR/(tCO$_2$/h)/km\]     |                    |                    |
+----------------------------+--------------------+--------------------+
| FOM \[%/yr\]               | 0.9                | 0.5                |
+----------------------------+--------------------+--------------------+
| Lifetime \[yr\]            | 50                 | 50                 |
+----------------------------+--------------------+--------------------+
| Network cost factor $f_y$  | 5.0 / 3.0 / 1.5 (2030 / 2040 / 2050)    |
+----------------------------+-----------------------------------------+

: CO~2~ pipeline cost parameters. Unit investment from the Danish Energy
Agency Technology Data for Energy Transport (2021), 12-inch pipeline
(120--500 tCO$_2$/h). The network cost factor $f_y$ scales total
pipeline capital cost by planning year.
:::

### CDR target {#subsubsec:app_policy_goals}

The policy goals set how much engineered removal (DACCS and BECCS)
European climate frameworks aim to incentivise in each scenario year.
They are not imposed as constraints in the optimisation. Instead, the
model sweeps the standalone CDR credit price and records the
cost-optimal deployment at each price, tracing a supply curve. For each
year the policy goal then selects the target-binding credit price: the
lowest credit price at which deployed CDR meets the target volume. The
goals thus represent demand-side ambition, and the target-binding credit
price is the market outcome that satisfies it.

The trajectory covers EU-31 (EU-27 plus Iceland, Liechtenstein, and
Norway, the NEGEM project scope) with the United Kingdom added
separately. Switzerland is omitted as negligible ($<0.5$ MtCO~2~/yr).
The EU-31 estimate is taken from the NEGEM expert elicitation under the
IEA Net-Zero Emissions scenario, which synthesises 34 expert interviews
into best estimates for each horizon [@david_reiner_quantifying_2023].
The UK contribution is the national total engineered greenhouse-gas
removal: the Government's 2030 target [@uk_government_net_2022] and the
Climate Change Committee's Seventh Carbon Budget balanced pathway for
2040 and 2050 [@bunting_analysis_2025].

::: {#tab:cdr-demand}
  **Year**                 **DACCS (EU-31)**               **BECCS (EU-31)**                  **EU-31 total**   **UK GGR**   **Adopted goal**
  ---------- ------------------------------- ------------------------------- -------------------------------- ------------ ------------------
  2030         $\sim$`<!-- -->`{=html}0.5--2    $\sim$`<!-- -->`{=html}5--15   $\sim$`<!-- -->`{=html}5.5--17            5                 22
  2040         $\sim$`<!-- -->`{=html}50--80   $\sim$`<!-- -->`{=html}40--60   $\sim$`<!-- -->`{=html}90--140           20                160
  2050                                   353                             131                              484           30                514

  : Engineered-CDR deployment goals and the values adopted as model
  policy goals. The adopted goal is the upper bound of the EU-31
  estimate plus the UK total engineered GGR. Units: MtCO~2~/yr.
:::

The adopted goals use the upper bound of the EU-31 range, consistent
with the NZE scenario as an upper envelope of ambition rather than a
central expectation: 22 Mt in 2030 (17 plus 5), 160 Mt in 2040 (140 plus
20), and 514 Mt in 2050 (484 plus 30). DACCS estimates in particular
inherit wide expert uncertainty in 2030 and 2040. These goals bound
credited CDR only and are independent of the cost sensitivtiies in
Appendix [7.1.6](#subsubsec:app_tech_cost_scenarios){reference-type="ref"
reference="subsubsec:app_tech_cost_scenarios"}.

### EU ETS price trajectory {#subsubsec:app_ets_prices}

The main analysis fixes the EU ETS price on a single medium trajectory,
119/279/463 €/t for 2030/2040/2050. It follows the PRIMES EU ETS 1
projection of [@gunther_carbon_2025], held fixed across that study's
scenarios, and is adopted for its consistency with European policy
scenarios and its integration of energy-system dynamics. Its 2050 value
matches the independent CAKE projection of 460 €/t
[@boratynski_linking_2024]. The standalone credit cancels the prevailing
ETS carbon value on every eligible capture link
(Section [3.1.1](#subsubsec:supply-modules){reference-type="ref"
reference="subsubsec:supply-modules"}), so the credit price that clears
the CDR target should be largely independent of the ETS level. The low
and high trajectories that test this, with their sources, are given with
the robustness check in
Appendix [7.2.4](#app:sensitivity_ets){reference-type="ref"
reference="app:sensitivity_ets"}.

**Scenario 2: Financing Gap**\

------------------------------------------------------------------------

Scenario 2 compares the target-binding credit price from Scenario 1
against what buyers will pay. The two benchmarks below, stated and
revealed, give the willingness-to-pay bands used in the gap assessment.

### Willingness to pay: construction of the benchmarks {#app:wtp-construction}

This subsection details the anchor derivation, currency conversion,
projection rationale, and cross-checks behind the benchmarks summarised
in Section [3.1.3](#subsubsec:wtp){reference-type="ref"
reference="subsubsec:wtp"} and reported in
Table [5](#tab:wtp-benchmarks){reference-type="ref"
reference="tab:wtp-benchmarks"}, with the full USD and EUR breakdown
across the three horizons in
Table [18](#tab:wtp-usd-eur){reference-type="ref"
reference="tab:wtp-usd-eur"}. The benchmarks draw on the 2026
CDR.fyi/OPIS buyer survey for the stated path and on Frontier and
CDR.fyi transaction records for the revealed path, cross-checked against
the EU ETS allowance level and Sweden's Klimatpremie clearing prices.
The uniform figure in each is the equal-weight mean of the BECCS and
DACCS values, so DACCS, with its higher unit cost and lower near-term
volume, is not diluted in a volume-weighted mean.

#### Rationale for two benchmarks.

The two benchmarks are imperfect proxies for the same quantity, the
price buyers will pay for durable CDR, with opposite and complementary
weaknesses. The survey is broad (all buyers, both methods, a forward
2030 value) but elicits intentions, which are documented to exceed
realised prices [@rodemeier_willingness_2023]. Transactions record
committed prices but are sparse, concentrated among a few early buyers,
and carry a first-of-a-kind premium. Reporting both and comparing each
against the modelled credit price separately tests whether the financing
gap closes under a survey-based or a transaction-based signal rather
than resting on either alone.

#### Stated anchors.

The 2030 medium is the arithmetic mean of the survey's method-level
buyer values for BECCS and DACCS, \$217 and \$361, giving \$289
(Table [16](#tab:survey-methods){reference-type="ref"
reference="tab:survey-methods"}). The high and low sensitivities apply a
proportional margin of $\pm 21\%$, giving \$350 and \$228. The high
coincides with the blended supplier reasonable-profit price for the two
methods ($\approx$\$340), an independent sell-side cross-check. A
proportional margin is used because buyer price tolerance scales with
the price level and is more stable across survey editions than the
absolute level. The $21\%$ width is calibrated from the spread between
the "good value" (\$208) and "expensive" (\$320) valuations elicited in
the 2025 edition (the only edition asking for those bounds, midpoint
\$264) and applied to the 2026 medium. At $0.864$ EUR/USD the 2030
stated anchors are €197.0, €249.7, and €302.4.

::: {#tab:survey-methods}
  Method     Buyer WTP   Supplier reasonable profit   Gap
  -------- ----------- ---------------------------- -----
  BECCS            217                          279    62
  DACCS            361                          400    39

  : 2026 CDR.fyi/OPIS survey: 2030 buyer WTP and supplier
  reasonable-profit prices for the two methods used in the benchmark
  (USD/tCO~2~). The supplier values give the $\approx$\$340 high-anchor
  cross-check.
:::

Source: CDR.fyi and OPIS (2026), 2030 expectations.

#### Revealed anchors.

Each method value is read from the deal-level register below
(Tables [\[tab:derived\]](#tab:derived){reference-type="ref"
reference="tab:derived"}
and [\[tab:disclosed\]](#tab:disclosed){reference-type="ref"
reference="tab:disclosed"}), and the two methods are then blended with
equal weight. The medium is the blend of the per-method means, \$316
(BECCS) and \$566 (DACCS), giving \$441. The low is the blend of the
per-method minima, \$214 (BECCS, the CO280 pulp-and-paper offtake) and
\$100 (DACCS, the Google--Holocene deal net of the United States 45Q tax
credit), giving \$157. The high is the blend of the per-method maxima,
\$427 (BECCS, Reverion) and \$989 (DACCS, Heirloom), giving \$708, with
a single \$2,345 DACCS order excluded as unrepresentative. At
$0.864$ EUR/USD the 2030 revealed anchors are €136, €381, and €612.

These revealed anchors rest on a deal-level register, assembled for
reproducibility. Every per-tonne price is either derived by formula from
a public contract value and tonnage (Block A) or taken from a direct
public per-tonne disclosure (Block B). Contracts whose price is only
estimated or undisclosed are excluded from the priced figures.

#### Derivation formula (Block A).

For every derived deal the implied unit price is obtained as
$$\begin{equation}
\text{Price}\ (\$/tCO\textsubscript{2}) \;=\; \frac{\text{Contract value (\$)} \times 10^{6}}
{\text{Contracted volume (tCO\textsubscript{2}{})}},
\label{eq:derive}
\end{equation}$$ where contract value is expressed in millions of US
dollars. The calculation is implemented as a spreadsheet formula rather
than a hand-entered number, so the derived column recomputes
automatically if an input is corrected during verification.

#### Validation.

Two independent checks support the derived figures. First, the values
produced by Equation [\[eq:derive\]](#eq:derive){reference-type="ref"
reference="eq:derive"} reproduce the \"calc price\" column of the
author's advisory deck deal-for-deal, cross-checking both inputs and
arithmetic. Second, each Frontier deal traces to a Frontier
"Writing" / portfolio page and, where available, to a CDR.fyi market
recap, so the contract value and tonnage entering the formula are
themselves primary-sourced.

#### Block A: derived transaction prices.

Block A comprises eighteen Frontier offtakes (F01--F18) and the
Google--Holocene bilateral deal (F19), with the per-tonne price computed
via Equation [\[eq:derive\]](#eq:derive){reference-type="ref"
reference="eq:derive"}. Values are in USD
(Table [\[tab:derived\]](#tab:derived){reference-type="ref"
reference="tab:derived"}).

::: longtable
\@L1.0cm L2.9cm L3.0cm C1.4cm C1.8cm C1.3cm@

\
**ID** & **Supplier** & **Method** & **Value (\$M)** & **Volume
(tCO~2~)** & **\$/tCO~2~**\
\
**ID** & **Supplier** & **Method** & **Value (\$M)** & **Volume
(tCO~2~)** & **\$/tCO~2~**\
\
F01 & NULIFE & BiCRS (bio-oil) & 44.2 & 122,000 & 362\
F02 & Reverion & BECCS (biogas) & 41.0 & 96,000 & 427\
F03 & Planetary & Ocean alkalinity & 31.3 & 115,211 & 272\
F04 & Arbor & BECCS & 41.0 & 116,000 & 353\
F05 & Hafslund Celsio & BECCS (WtE) & 31.6 & 100,000 & 316\
F06 & Eion & Enhanced weathering & 33.0 & 78,707 & 419\
F07 & Phlair & DAC (electrochem.) & 30.6 & 47,000 & 651\
F08 & CO280 & BECCS (pulp & paper) & 48.0 & 224,500 & 214\
F09 & CREW Carbon & Alkalinity (wastewater) & 32.1 & 71,878 & 447\
F10 & Terradot & Enhanced weathering & 27.0 & 90,000 & 300\
F11 & CarbonRun & River alkalinity & 25.4 & 55,442 & 458\
F12 & 280 Earth & DAC & 40.0 & 61,571 & 650\
F13 & Stockholm Exergi & BECCS & 48.6 & 179,998 & 270\
F14 & Vaulted Deep & Biomass storage & 58.3 & 152,000 & 384\
F15 & Heirloom & DAC & 26.6 & 26,889 & 989\
F16 & CC & Heirloom & DAC & 20.0 & 45,500 & 440\
F17 & Charm Industrial & BiCRS (bio-oil) & 53.0 & 112,000 & 473\
F18 & Lithos Carbon & Enhanced weathering & 57.1 & 154,240 & 370\
F19 & Holocene & DAC (Google) & 10.0 & 100,000 & 100^a^\
:::

^a^ Lowest disclosed DAC price. It is net of the US 45Q tax credit
(\~\$180/tCO~2~) on an unsubsidised \~\$280/tCO~2~, and is a forward
learning-curve sale. Not treated as a market or deployment cost.

#### Unweighted method averages.

Simple (unweighted) means of the Block A prices give indicative
reference points by method: DAC \$566/tCO~2~, BECCS \$316/tCO~2~, and
enhanced weathering \$363/tCO~2~. These are means of *forward* offtake
contracts and therefore sit above current spot levels. They are reported
with that caveat and are not deployment-cost estimates. The DAC mean
includes the subsidy-net Holocene point.

#### Block B: disclosed per-tonne prices.

Block B records deals for which a buyer disclosed a per-tonne price
directly, without a value/tonnage split, so
Equation [\[eq:derive\]](#eq:derive){reference-type="ref"
reference="eq:derive"} does not apply
(Table [\[tab:disclosed\]](#tab:disclosed){reference-type="ref"
reference="tab:disclosed"}). The figure is taken from the public
disclosure and cross-checked against CDR.fyi market recaps. The
CO280--JPMorgan contract is one of the lowest disclosed engineered-CDR
prices: biogenic capture at a US Gulf Coast pulp-and-paper mill,
contracted at under \$200/tCO~2~. It is distinct from the
Frontier--CO280 offtake (F08, 224,500 tCO~2~) and from the larger
Microsoft--CO280 contract, which are recorded separately.

::: tabular
L1.0cm L1.8cm L2.6cm L2.4cm C1.8cm C1.6cm **ID** & **Supplier** &
**Method** & **Buyer** & **Volume (tCO~2~)** & **\$/tCO~2~**\
U01 & CO280 & BECCS (pulp & paper) & JPMorgan Chase & 450,000 &
$<\!200$\
:::

#### Currency conversion.

A fixed rate of $1$ USD${}={}$$0.864$ EUR (EUR/USD $\approx 1.157$) is
held constant across years, so the conversion adds no exchange-rate
movement to the trajectory. Underlying data and ratios remain in USD. A
$\pm 5\%$ variation in the rate shifts the 2030 anchors by roughly €10
to €16 per tonne, small relative to the spread between sensitivities.

#### Projection to 2040 and 2050.

Neither source reports beyond 2030, so the anchors are projected at
scenario-specific annual rates. For the stated benchmark: low
$-3.5\%$/yr (learning and commoditisation), medium $-1\%$/yr (cost
softening partly offset by rising carbon prices), high $+2.5\%$/yr
(compliance value and scarcity), the last rising toward but below the
modelled EU ETS price
(Appendix [7.1.12](#subsubsec:app_ets_prices){reference-type="ref"
reference="subsubsec:app_ets_prices"}). For the revealed benchmark:
medium $\approx-1.1\%$/yr and low $\approx-0.9\%$/yr, reflecting the
unwinding of the first-of-a-kind premium and convergence toward the
stated level (survey buyer-supplier gap \$107 to \$98 to \$48). The
revealed high rises at $+2.5\%$/yr on the premise that durable removals
are not fungible with allowances and may price above them
[@rickels_integrating_2021]. The downward convergence is also consistent
with the demand curve: the survey ($n=20$ buyers) and today's deals sit
at its early, high-value end, so reaching the EU deployment goals
(Appendix [7.1.11](#subsubsec:app_policy_goals){reference-type="ref"
reference="subsubsec:app_policy_goals"}) draws in lower-value marginal
buyers and lowers the clearing price
[@yang_global_2023; @johnstone_carbon_2026]. Falling technology costs
reinforce the medium and low paths
[@danish_energy_agency_technology_2024].

#### Cross-checks.

The constructed ranges are consistent with two independent sources.
First, observed BECCS offtake prices fall sharply with contract duration
(Table [17](#tab:beccs-duration){reference-type="ref"
reference="tab:beccs-duration"}): long-duration committed deals (10--15
years), typical of capital-intensive durable CDR, cluster at roughly
\$100--\$130/t (€92--€120), sitting just below the revealed low anchor
(€136) and the documented discount of committed buyers relative to
survey-stated WTP [@frontier_climate_frontier_2025]. Second, two
compliance benchmarks bracket the stated range. The near-term EU ETS
allowance level ($\approx$€60--€80/t) sits below all three 2030 anchors,
confirming that buyers would prefer removals to allowances in the near
term [@rickels_integrating_2021]. Clearing prices from Sweden's
Klimatpremie reverse auction ($\approx$€100--€270/t) bracket the stated
medium (€249), the only operational European benchmark for
government-backed CDR procurement
[@lundberg_missing_2022; @fridahl_potential_2024].

::: {#tab:beccs-duration}
  Contract duration     Approximate price (USD/tCO~2~)
  ------------------- --------------------------------
  3 years                                     260--353
  5 years                                          195
  10 years                                         115
  15 years                                         100

  : Observed BECCS offtake prices by contract duration, used to
  cross-check the revealed anchors.
:::

Source: Frontier Climate (2025) and CDR.fyi transaction data. Prices
reflect publicly disclosed offtakes. Feedstock, capture technology, and
storage location also influence pricing.

#### Limitations.

First, the stated benchmark rests on a small buyer sample (20 valid
purchaser responses) and reflects intentions rather than executed
behaviour. Stated valuations exceed revealed in level even where their
direction is informative [@rodemeier_willingness_2023], consistent with
the revealed benchmark exceeding the stated medium in 2030 only because
of the early-deal premium, not a level transfer. Second, the 2040 and
2050 values are projections under assumed rates, reported as a
sensitivity range, not forecasts. Third, the uniform price aggregates
two pathways with different cost structures, so the method-level values
are retained in the source register.

::: {#tab:wtp-usd-eur}
+-------------+-----------------------------+-----------------------------+
|             | USD                         | EUR                         |
+:============+========:+========:+========:+========:+========:+========:+
| 2-4 (lr)5-7 | 2030    | 2040    | 2050    | 2030    | 2040    | 2050    |
| Sensitivity |         |         |         |         |         |         |
+-------------+---------+---------+---------+---------+---------+---------+
| *Stated preference*                                                     |
+-------------+---------+---------+---------+---------+---------+---------+
| high        | 350     | 448     | 573     | 302.4   | 387.1   | 495.1   |
+-------------+---------+---------+---------+---------+---------+---------+
| medium      | 289     | 261     | 236     | 249.7   | 225.5   | 203.9   |
+-------------+---------+---------+---------+---------+---------+---------+
| low         | 228     | 160     | 112     | 197.0   | 138.2   | 96.8    |
+-------------+---------+---------+---------+---------+---------+---------+
|             |         |         |         |         |         |         |
+-------------+---------+---------+---------+---------+---------+---------+
| high        | 708     | 906     | 1,160   | 612     | 783     | 1,002   |
+-------------+---------+---------+---------+---------+---------+---------+
| medium      | 441     | 395     | 353     | 381     | 341     | 305     |
+-------------+---------+---------+---------+---------+---------+---------+
| low         | 157     | 143     | 131     | 136     | 124     | 113     |
+-------------+---------+---------+---------+---------+---------+---------+

: Willingness-to-pay benchmarks in USD and EUR (per tonne CO~2~).
:::

EUR figures rounded. Minor differences from direct conversion reflect
rounding.

## Results {#subsec:app_results}

### Baseline {#subsec:app_results_baseline}

Figure [13](#fig:baseline_generation_mix){reference-type="ref"
reference="fig:baseline_generation_mix"} shows the baseline electricity
mix.

<figure id="fig:baseline_generation_mix" data-latex-placement="H">
<img src="./Figures/Appendix/baseline_energymix.png" />
<figcaption>Baseline electricity generation mix in the medium-cost
zero-credit case.</figcaption>
</figure>

Total generation rises from 3962 TWh/yr in 2030 to 6194 TWh/yr in 2050
as electrification spreads across heat, transport, and industry. Wind
and solar carry that growth. Wind rises from 36% to 53% of supply and
solar from 24% to 36%, so the two variable renewables together move from
60% to 89% over the period. Hydro holds near 471 TWh/yr throughout, a
fixed resource that shrinks in share as demand grows. Gas is small and
falls to 0.8% by 2050, running only as the flexible remainder priced by
the ETS, since no emissions cap is applied.

Nuclear falls from 952 TWh/yr (24%) in 2030 to 213 TWh/yr (4%) in 2040
and 133 TWh/yr (2%) in 2050. The decline is economic, not a policy
constraint. The country ceilings in
Table [8](#tab:app_nuclear_limits){reference-type="ref"
reference="tab:app_nuclear_limits"} are upper bounds with no minimum, so
the optimiser builds nuclear only where it beats the alternative. In
2030 the existing fleet still runs and supplies a quarter of demand. As
that fleet retires on its lifetime schedule, the model does not rebuild
it, because wind and solar deliver the same energy at lower cost. The
ceilings even rise over time, with France reaching 82 GW by 2050, but
the cost-optimal system leaves most of that headroom unused. The Belgian
and Spanish phase-outs remove part of the fleet, yet the main driver is
non-replacement of retiring capacity rather than mandated closure.

#### Spatial pattern of baseline capture.

Tables [19](#tab:app_baseline_spatial_2030){reference-type="ref"
reference="tab:app_baseline_spatial_2030"}
to [21](#tab:app_baseline_spatial_2050){reference-type="ref"
reference="tab:app_baseline_spatial_2050"} give the country detail
behind the baseline map
(Figure [3](#fig:baseline_map){reference-type="ref"
reference="fig:baseline_map"}), for countries above 1 Mt /yr of gross
capture. All baseline capture is fossil or process . No BECCS appears at
zero credit price, and the only CDR is 11.3 Mt /yr of DACCS in 2050.
Captured is gross routed off the point source, Stored is geological
sequestration, and Utilised is sent to synthetic fuels and later
re-emitted. In 2030 and 2040 capture and storage track each other
closely and stay near the capture site. By 2050 most captured goes to
synthetic fuels instead of storage, so the Utilised column carries most
of the volume and sequestration falls to 24 Mt /yr.

::: {#tab:app_baseline_spatial_2030}
  Country                 Captured   Stored
  --------------------- ---------- --------
  United Kingdom (GB)         22.4     22.4
  Norway (NO)                 11.3     11.3
  Netherlands (NL)             9.1      9.1
  Denmark (DK)                 2.5      2.5
  Germany (DE)                 1.3      1.3
  EU total                      47       47

  : Baseline per-country capture and fate, medium cost, 2030. All
  capture is fossil or process . Mt /yr.
:::

::: {#tab:app_baseline_spatial_2040}
  Country                 Captured   Stored   Utilised
  --------------------- ---------- -------- ----------
  United Kingdom (GB)         16.7     16.2        0.5
  Netherlands (NL)            13.2     14.1        0.0
  Norway (NO)                 10.8     10.8        0.0
  Italy (IT)                  10.1      8.7        1.3
  Spain (ES)                   8.7      2.9        5.7
  Germany (DE)                 7.6      7.6        0.0
  Portugal (PT)                4.4      1.3        3.1
  Ireland (IE)                 4.4      1.6        2.8
  Denmark (DK)                 3.0      3.0        0.0
  Lithuania (LT)               2.1      2.1        0.0
  Sweden (SE)                  2.0      2.0        0.0
  Poland (PL)                  1.5      1.5        0.0
  EU total                      87       73         14

  : Baseline per-country capture and fate, medium cost, 2040. All
  capture is fossil or process . Mt /yr.
:::

::: {#tab:app_baseline_spatial_2050}
  Country                 Fossil & process   DACCS   Captured   Stored   Utilised
  --------------------- ------------------ ------- ---------- -------- ----------
  Spain (ES)                          19.4     8.9       28.3      2.2       26.1
  Italy (IT)                          13.9     1.5       15.4      3.3       12.1
  United Kingdom (GB)                 11.6     0.1       11.7      2.1        9.7
  Norway (NO)                         10.7     0.0       10.7      1.9        8.8
  Germany (DE)                        10.2     0.0       10.2      2.5        8.0
  Netherlands (NL)                     8.8     0.0        8.8      5.4        2.2
  France (FR)                          4.9     0.2        5.1      0.0        5.4
  Greece (GR)                          4.5     0.0        4.5      1.0        3.5
  Portugal (PT)                        4.4     0.0        4.4      0.9        3.5
  Ireland (IE)                         4.2     0.0        4.2      0.7        3.4
  Poland (PL)                          4.0     0.0        4.0      1.5        2.4
  EU total                             116      11        127       24        103

  : Baseline per-country capture and fate, medium cost, 2050. DACCS is
  the only CDR. Mt /yr.
:::

### Scenario 1: deployment response {#app:scenario_tables}

#### Credit-price sweep

<figure id="app:creditpricesweep" data-latex-placement="H">
<img src="./Figures/Appendix/fig_credit_price_sweep_medium.png" />
<figcaption>Credit price sweep in the medium cost
sensitivitiy</figcaption>
</figure>

#### Capture mix and fate.

The maps in Figures [3](#fig:baseline_map){reference-type="ref"
reference="fig:baseline_map"} and [5](#fig:cdr_map){reference-type="ref"
reference="fig:cdr_map"} resolve capture and storage by country but
cannot show the full technology and routing split for every
scenario--year. Figure [15](#fig:app_capture_fate){reference-type="ref"
reference="fig:app_capture_fate"} gives the EU-aggregate companion,
computed from the same nodal energy balances and sequestration-link
flows. *Baseline* is the medium-cost, zero-credit case. *low*, *medium*,
and *high* are each taken at that cost sensitivity's target-binding
clearing price. Shading runs from pale (low share) to deep red (high
share) on a common 0--100 % scale: the left panel shows the entry order
of engineered removal (BECCS first, DAC overtaking from 2040), the right
panel the shift from utilisation to permanent storage.

<figure id="fig:app_capture_fate" data-latex-placement="H">
<div class="minipage">
<p><span>Capture mix (% of gross capture)</span><br />
</p>
<table>
<tbody>
<tr>
<td style="text-align: left;">Year</td>
<td style="text-align: left;">Scenario</td>
<td style="text-align: right;">BECCS</td>
<td style="text-align: right;">DAC</td>
<td style="text-align: right;">Ind. CCS</td>
<td style="text-align: right;">Capture</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: right;">[%]</td>
<td style="text-align: right;">[%]</td>
<td style="text-align: right;">[%]</td>
<td style="text-align: right;">[Mt/yr]</td>
</tr>
<tr>
<td style="text-align: left;">2030</td>
<td style="text-align: left;">Baseline</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;"><span
style="color: white">100.0</span></td>
<td style="text-align: right;">47</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">low</td>
<td style="text-align: right;">35.1</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;"><span
style="color: white">64.9</span></td>
<td style="text-align: right;">70</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">medium</td>
<td style="text-align: right;">33.3</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;"><span
style="color: white">66.7</span></td>
<td style="text-align: right;">70</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">high</td>
<td style="text-align: right;"><span
style="color: white">69.8</span></td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">30.2</td>
<td style="text-align: right;">70</td>
</tr>
<tr>
<td style="text-align: left;"><span>1-6</span> 2040</td>
<td style="text-align: left;">Baseline</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;"><span
style="color: white">100.0</span></td>
<td style="text-align: right;">87</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">low</td>
<td style="text-align: right;">10.0</td>
<td style="text-align: right;"><span
style="color: white">64.8</span></td>
<td style="text-align: right;">25.2</td>
<td style="text-align: right;">335</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">medium</td>
<td style="text-align: right;">20.0</td>
<td style="text-align: right;">61.9</td>
<td style="text-align: right;">18.2</td>
<td style="text-align: right;">337</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">high</td>
<td style="text-align: right;">25.3</td>
<td style="text-align: right;">53.8</td>
<td style="text-align: right;">21.0</td>
<td style="text-align: right;">337</td>
</tr>
<tr>
<td style="text-align: left;"><span>1-6</span> 2050</td>
<td style="text-align: left;">Baseline</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">8.9</td>
<td style="text-align: right;"><span
style="color: white">91.1</span></td>
<td style="text-align: right;">126</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">low</td>
<td style="text-align: right;">4.8</td>
<td style="text-align: right;"><span
style="color: white">79.3</span></td>
<td style="text-align: right;">15.9</td>
<td style="text-align: right;">704</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">medium</td>
<td style="text-align: right;">4.6</td>
<td style="text-align: right;"><span
style="color: white">78.8</span></td>
<td style="text-align: right;">16.5</td>
<td style="text-align: right;">699</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">high</td>
<td style="text-align: right;">12.3</td>
<td style="text-align: right;"><span
style="color: white">72.8</span></td>
<td style="text-align: right;">14.9</td>
<td style="text-align: right;">700</td>
</tr>
</tbody>
</table>
</div>
<div class="minipage">
<p><span> fate (% of routed )</span><br />
</p>
<table>
<tbody>
<tr>
<td style="text-align: left;">Year</td>
<td style="text-align: left;">Scenario</td>
<td style="text-align: right;">Sequest.</td>
<td style="text-align: right;">Util.</td>
<td style="text-align: right;">Routed</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: right;">[%]</td>
<td style="text-align: right;">[%]</td>
<td style="text-align: right;">[Mt/yr]</td>
</tr>
<tr>
<td style="text-align: left;">2030</td>
<td style="text-align: left;">Baseline</td>
<td style="text-align: right;"><span
style="color: white">100.0</span></td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">47</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">low</td>
<td style="text-align: right;"><span
style="color: white">100.0</span></td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">70</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">medium</td>
<td style="text-align: right;"><span
style="color: white">100.0</span></td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">70</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">high</td>
<td style="text-align: right;"><span
style="color: white">100.0</span></td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">70</td>
</tr>
<tr>
<td style="text-align: left;"><span>1-5</span> 2040</td>
<td style="text-align: left;">Baseline</td>
<td style="text-align: right;"><span
style="color: white">84.4</span></td>
<td style="text-align: right;">15.6</td>
<td style="text-align: right;">87</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">low</td>
<td style="text-align: right;"><span
style="color: white">95.5</span></td>
<td style="text-align: right;">4.5</td>
<td style="text-align: right;">335</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">medium</td>
<td style="text-align: right;"><span
style="color: white">95.0</span></td>
<td style="text-align: right;">5.0</td>
<td style="text-align: right;">337</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">high</td>
<td style="text-align: right;"><span
style="color: white">94.9</span></td>
<td style="text-align: right;">5.1</td>
<td style="text-align: right;">337</td>
</tr>
<tr>
<td style="text-align: left;"><span>1-5</span> 2050</td>
<td style="text-align: left;">Baseline</td>
<td style="text-align: right;">19.0</td>
<td style="text-align: right;"><span
style="color: white">81.0</span></td>
<td style="text-align: right;">126</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">low</td>
<td style="text-align: right;"><span
style="color: white">85.3</span></td>
<td style="text-align: right;">14.7</td>
<td style="text-align: right;">704</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">medium</td>
<td style="text-align: right;"><span
style="color: white">85.8</span></td>
<td style="text-align: right;">14.2</td>
<td style="text-align: right;">699</td>
</tr>
<tr>
<td style="text-align: left;"></td>
<td style="text-align: left;">high</td>
<td style="text-align: right;"><span
style="color: white">85.7</span></td>
<td style="text-align: right;">14.3</td>
<td style="text-align: right;">700</td>
</tr>
</tbody>
</table>
</div>
<figcaption>EU-aggregate capture mix (left) and fate (right) at the
target-binding credit price, by horizon and cost sensitivity. Shares in
percent (rows sum to 100), final column absolute (Mt /yr). Baseline is
medium cost, zero credit.</figcaption>
</figure>

#### Headline metrics across cost sensitivities.

Table [22](#tab:app_cost_headline){reference-type="ref"
reference="tab:app_cost_headline"} shows target-binding credit price
across the low, medium, and high cost sensitivities, with one column
block per sensitivity over the three horizons. *Captured* is the gross
taken from each source (DACCS, BECCS, and fossil or process CCS); this
captured either reaches geological storage (*Stored*) or is routed to
synthetic fuels (*Utilised*). *Credited CDR* is the narrower quantity
that earns a credit: the stored, non-fossil share, that is, the DAC and
biogenic that is permanently sequestered. It therefore excludes both
utilised and all fossil & process capture.

&

::: {#tab:app_cost_headline}
+--------------------+--------------------+--------------------+--------------------+
|                    | low                | medium             | high               |
+:===================+:====:+:====:+:====:+:====:+:====:+:====:+:====:+:====:+:====:+
| 2-4(lr)5-7(lr)8-10 | 2030 | 2040 | 2050 | 2030 | 2040 | 2050 | 2030 | 2040 | 2050 |
| Metric             |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Binding credit     | 100  | 100  | 100  | 100  | 150  | 100  | 150  | 150  | 150  |
| price (€/t)        |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Total Captured     | 70   | 335  | 704  | 70   | 336  | 699  | 70   | 337  | 700  |
| (Mt /yr)           |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Captured DACCS     | 0    | 217  | 558  | 0    | 208  | 551  | 0    | 181  | 509  |
| (Mt /yr)           |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Captured BECCS     | 25   | 34   | 34   | 23   | 67   | 33   | 49   | 85   | 86   |
| (Mt /yr)           |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Fossil & process   | 45   | 84   | 112  | 47   | 61   | 115  | 21   | 71   | 105  |
| CCS (Mt /yr)       |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Stored (Mt /yr)    | 70   | 320  | 600  | 70   | 320  | 600  | 70   | 320  | 600  |
+--------------------+------+------+------+------+------+------+------+------+------+
| Utilised (Mt /yr)  | 0    | 15   | 104  | 0    | 17   | 99   | 0    | 17   | 100  |
+--------------------+------+------+------+------+------+------+------+------+------+
| Credited CDR       | 25   | 251  | 592  | 23   | 275  | 584  | 49   | 266  | 595  |
| (Mt /yr)           |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| of which DACCS     | 0    | 217  | 558  | 0    | 208  | 551  | 0    | 181  | 509  |
+--------------------+------+------+------+------+------+------+------+------+------+
| of which BECCS     | 25   | 34   | 34   | 23   | 67   | 33   | 49   | 85   | 86   |
+--------------------+------+------+------+------+------+------+------+------+------+
| Average            | 68   | 99   | 67   | 68   | 98   | 69   | 68   | 100  | 68   |
| electricity price  |      |      |      |      |      |      |      |      |      |
| (€/MWh)            |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Average heat price | 40   | 91   | 48   | 40   | 92   | 52   | 40   | 94   | 50   |
| (€/MWh)            |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| $\Delta$ system    | 0    | 15   | 33   | 2    | 23   | 46   | 3    | 32   | 66   |
| cost (bn €/yr)     |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+

: Headline metrics at the target-binding credit price by cost
sensitivity, each block one sensitivity over 2030/2040/2050. Captured
rows are gross capture by origin; captured is either geologically stored
or utilised (routed to synthetic fuels). Credited CDR is the stored,
non-fossil (DAC + BECCS) share of capture under waterfall attribution,
split by origin.
:::

#### Spatial distribution of hosting.

Tables [24](#tab:app_s1_spatial_2040){reference-type="ref"
reference="tab:app_s1_spatial_2040"}
and [25](#tab:app_s1_spatial_2050){reference-type="ref"
reference="tab:app_s1_spatial_2050"} give the per-country capture behind
the scenario 1 map (Figure [5](#fig:cdr_map){reference-type="ref"
reference="fig:cdr_map"}), read at the target-binding credit price
(100 €/t in 2030, 150 €/t in 2040 and 2050) in the medium-cost case. In
2030 all credited CDR is BECCS and stays small, so only the later years
are tabulated. The columns split gross capture into BECCS, DACCS, and
competing fossil or process CCS, next to geological storage and
utilisation. The Captured and Stored columns match the hosting columns
of Tables [31](#tab:app_country_cost_2040){reference-type="ref"
reference="tab:app_country_cost_2040"}
and [32](#tab:app_country_cost_2050){reference-type="ref"
reference="tab:app_country_cost_2050"}. Capture concentrates in Denmark
and the United Kingdom, which together hold 65% of capture in 2040 and
68% in 2050. DACCS drives that concentration, while BECCS stays spread
across biomass-rich countries such as Sweden, Italy, and Portugal. Spain
and Italy capture more than they store and route the surplus to
synthetic fuels, which is why their Utilised column is large.

::: {#tab:app_s1_spatial_2030}
  Country                 BECCS   Fossil & proc.   Captured   Stored
  --------------------- ------- ---------------- ---------- --------
  United Kingdom (GB)       7.2             22.4       29.7     29.7
  Norway (NO)               1.7             11.3       13.0     13.0
  Netherlands (NL)          0.1              8.9        9.1      9.1
  Sweden (SE)               5.0              0.0        5.0      5.0
  Denmark (DK)              1.6              2.7        4.3      4.3
  Spain (ES)                2.9              0.0        2.9      2.9
  Italy (IT)                2.7              0.0        2.7      2.7
  Germany (DE)              1.0              1.3        2.3      2.3
  Portugal (PT)             1.0              0.0        1.0      1.0
  EU total                   23               47         70       70

  : Scenario 1 per-country capture and fate at the target-binding credit
  price (100 €/t), medium cost, 2030. No DACCS or utilisation appears
  this year. Mt /yr.
:::

::: {#tab:app_s1_spatial_2040}
  Country                 BECCS   DACCS   Fossil & proc.   Captured   Stored   Utilised
  --------------------- ------- ------- ---------------- ---------- -------- ----------
  Denmark (DK)              2.1   147.5              2.1      151.7    159.6        0.0
  United Kingdom (GB)       9.5    41.9             15.6       67.0     66.4        0.7
  Spain (ES)                6.6    10.2              2.9       19.7      5.9       13.8
  Sweden (SE)              16.1     0.0              0.5       16.5      8.7        0.0
  Germany (DE)              6.5     0.0              9.6       16.1     16.6        0.0
  Ireland (IE)              1.3     8.5              4.4       14.2     14.1        0.1
  Norway (NO)               2.3     0.0             10.1       12.4     12.4        0.0
  Netherlands (NL)          2.5     0.0              9.1       11.6     14.5        0.0
  Italy (IT)                5.5     0.1              4.4       10.0      7.9        2.1
  Portugal (PT)             5.3     0.0              0.0        5.3      5.3        0.0
  Poland (PL)               1.9     0.0              1.1        3.0      3.0        0.0
  EU total                   67     208               62        337      320         17

  : Scenario 1 per-country capture and fate at the target-binding credit
  price, medium cost, 2040. Mt /yr.
:::

::: {#tab:app_s1_spatial_2050}
  Country                 BECCS   DACCS   Fossil & proc.   Captured   Stored   Utilised
  --------------------- ------- ------- ---------------- ---------- -------- ----------
  United Kingdom (GB)       7.2   270.9             10.8      288.9    282.5        6.4
  Denmark (DK)              2.1   193.1              2.3      197.5    197.8        0.0
  Spain (ES)                3.8    35.7             10.8       50.4      8.2       42.2
  Italy (IT)                5.7    12.4             18.0       36.1     15.0       21.0
  Norway (NO)               2.3    15.7              9.5       27.5     27.2        0.4
  Ireland (IE)              0.1    14.3              4.2       18.7     10.7        8.0
  Germany (DE)              1.7     0.0             16.6       18.3     16.3        0.0
  Sweden (SE)               8.3     0.4              2.1       10.8     10.3        0.1
  Netherlands (NL)          0.1     0.0              9.7        9.9     18.3        0.0
  Greece (GR)               0.0     2.6              4.5        7.1      0.9        6.2
  France (FR)               0.0     3.0              3.4        6.4      0.0        5.8
  Belgium (BE)              0.0     0.0              5.8        5.8      0.0        0.0
  Portugal (PT)             1.0     0.0              4.4        5.5      5.5        0.0
  EU total                   32     551              116        699      600         99

  : Scenario 1 per-country capture and fate at the target-binding credit
  price, medium cost, 2050. Mt /yr.
:::

#### Robustness check of DAC deployment to weather {#app:sensitivity_weather}

Module 2 scales S-DAC and L-DAC energy demand by country climate
(Appendix [7.1.7](#app:dac-weather){reference-type="ref"
reference="app:dac-weather"}). This check asks whether that scaling
changes the cost or the siting of removal. It changes neither
materially.

Weather lowers DAC energy rather than raising it. The
deployment-weighted S-DAC energy-cost factor is 0.87 in 2050, against
0.90 unweighted and 1.0 under the uniform assumption. Siting at scale
therefore draws favourable climates about 13% below the reference. The
best-to-worst spread across European climates is small. At medium energy
prices the S-DAC energy cost runs from $-6$ €/t in Norway to $+2$ €/t in
Greece, a spread of 8 €/t. The spread is 5 €/t at low and 11 €/t at high
energy prices. All three are below one 50 €/t price step, so weather
does not move the target-binding price.

Siting is set by storage, not climate. Across countries the rank
correlation between the S-DAC energy-cost factor and S-DAC capture is
near zero (Spearman $\rho = 0.02$, $p = 0.92$ in 2050). The correlation
with geological sequestration is strong instead ($\rho = 0.57$ to
$0.84$, $p < 0.001$). The hosts, Denmark, the United Kingdom, and
Norway, are storage-rich, not climate-favoured.

#### Country-level electricity and heat price effects {#app:host_price_effects}

The price effects summarised in the Results vary by country.
Tables [26](#tab:app_host_elec_price_delta){reference-type="ref"
reference="tab:app_host_elec_price_delta"}
and [27](#tab:app_host_heat_price_delta){reference-type="ref"
reference="tab:app_host_heat_price_delta"} give the electricity and heat
price change for the four largest hosts, relative to the zero-credit
baseline, at each year's target-binding credit price (100 €/t in 2030
and 2050, and 150 €/t in 2040). The change stays within a few euro per
megawatt-hour in most country-years. The larger values fall in 2040 and
2050, once DAC enters, and they are negative as often as positive. The
EU average row reproduces the system figures reported in the Results.

::: {#tab:app_host_elec_price_delta}
  Country                           2030      2040      2050
  ---------------------------- --------- --------- ---------
  Denmark (DK)                       0.0   $+11.9$   $-27.4$
  United Kingdom (GB)                0.0       0.0    $+0.1$
  Spain (ES)                     $+14.1$    $+8.4$   $+12.4$
  Norway (NO)$^{\mathrm{a}}$      $+0.9$   $-40.7$    $+4.0$
  EU average                      $+0.9$    $+1.8$    $+5.2$

  : Change in the average country electricity price relative to the
  zero-credit baseline, medium cost, at each year's target-binding
  credit price. Values are the unweighted mean nodal price, binding
  minus baseline, in €/MWh, the same metric as the system average in the
  Results.
:::

$^{\mathrm{a}}$In Norway in 2040 the average is set by a change in the
number of low- and negative-priced hours, not by a shift in the typical
price. The median nodal-price change is near zero.

::: {#tab:app_host_heat_price_delta}
  Country                   2030     2040     2050
  --------------------- -------- -------- --------
  Denmark (DK)            $+8.9$   $+1.6$   $+3.4$
  United Kingdom (GB)     $+0.2$   $+8.4$   $-4.2$
  Spain (ES)              $+0.9$   $+1.9$   $+1.1$
  Norway (NO)             $+0.4$   $+1.2$   $-4.1$
  EU average              $+0.3$   $+4.5$   $-2.2$

  : Change in the average country heat price relative to the zero-credit
  baseline, medium cost, at each year's target-binding credit price.
  Values are the load-weighted mean heat-bus price, binding minus
  baseline, in €/MWh.
:::

### Scenario 2: financing gap {#app:financing_gap}

#### Price gap across sensitivities.

Tables [28](#tab:app_gap_stated){reference-type="ref"
reference="tab:app_gap_stated"}
and [29](#tab:app_gap_revealed){reference-type="ref"
reference="tab:app_gap_revealed"} give the gap between the
target-binding credit price and willingness to pay for every cost and
WTP combination, by year. Positive is a financing gap, negative means
WTP covers the price. In 2030 the gap is essentially closed, with a
small $+14$ €/t showing on the revealed benchmark only under high cost
and low WTP. It opens under low WTP from 2040 and is widest in 2050
under low WTP at high cost ($+54$ €/t stated, $+37$ €/t revealed).

::: {#tab:app_gap_stated}
+----------+--------+-------------------------+
|          |        | WTP level               |
+:=========+:=======+======:+=======:+=======:+
| 3-5 Year | Cost   | low   | medium | high   |
+----------+--------+-------+--------+--------+
| 2030     | low    | $-96$ | $-149$ | $-202$ |
+----------+--------+-------+--------+--------+
|          | medium | $-96$ | $-149$ | $-202$ |
+----------+--------+-------+--------+--------+
|          | high   | $-46$ | $-99$  | $-152$ |
+----------+--------+-------+--------+--------+
| 2040     | low    | $-38$ | $-125$ | $-286$ |
+----------+--------+-------+--------+--------+
|          | medium | $+12$ | $-75$  | $-236$ |
+----------+--------+-------+--------+--------+
|          | high   | $+12$ | $-75$  | $-236$ |
+----------+--------+-------+--------+--------+
| 2050     | low    | $+4$  | $-103$ | $-394$ |
+----------+--------+-------+--------+--------+
|          | medium | $+4$  | $-103$ | $-394$ |
+----------+--------+-------+--------+--------+
|          | high   | $+54$ | $-53$  | $-344$ |
+----------+--------+-------+--------+--------+

: Price gap (target-binding credit price minus WTP, €/t ) by cost
sensitivity and WTP level, stated benchmark. Positive is a financing
gap.
:::

::: {#tab:app_gap_revealed}
+----------+--------+-------------------------+
|          |        | WTP level               |
+:=========+:=======+======:+=======:+=======:+
| 3-5 Year | Cost   | low   | medium | high   |
+----------+--------+-------+--------+--------+
| 2030     | low    | $-36$ | $-281$ | $-512$ |
+----------+--------+-------+--------+--------+
|          | medium | $-36$ | $-281$ | $-512$ |
+----------+--------+-------+--------+--------+
|          | high   | $+14$ | $-231$ | $-462$ |
+----------+--------+-------+--------+--------+
| 2040     | low    | $-24$ | $-241$ | $-683$ |
+----------+--------+-------+--------+--------+
|          | medium | $+26$ | $-191$ | $-633$ |
+----------+--------+-------+--------+--------+
|          | high   | $+26$ | $-191$ | $-633$ |
+----------+--------+-------+--------+--------+
| 2050     | low    | $-13$ | $-205$ | $-902$ |
+----------+--------+-------+--------+--------+
|          | medium | $-13$ | $-205$ | $-902$ |
+----------+--------+-------+--------+--------+
|          | high   | $+37$ | $-155$ | $-852$ |
+----------+--------+-------+--------+--------+

: Price gap (€/t ) by cost sensitivity and WTP level, revealed
benchmark. Sign convention as in
Table [28](#tab:app_gap_stated){reference-type="ref"
reference="tab:app_gap_stated"}.
:::

#### Country-level distribution of the system cost.

The aggregate comparison in the Results sets the credit revenue against
the system-cost increase at the European level. That increase is not
spread evenly across countries.
Tables [30](#tab:app_country_cost_2030){reference-type="ref"
reference="tab:app_country_cost_2030"}
to [32](#tab:app_country_cost_2050){reference-type="ref"
reference="tab:app_country_cost_2050"} decompose it by country, one
table per planning year, next to the carbon removal each country hosts.
Each table is read at that year's target-binding credit price, 100 €/t
in 2030, 150 €/t in 2040, and 100 €/t in 2050, against the medium-cost
zero-credit baseline. The CDR targets met at these prices are 22, 160,
and 514 Mt /yr. The tables list countries that carry a net cost increase
above 0.05 bn €/yr. Countries below this reporting threshold are omitted
and may show small positive or negative changes.

The Total column is the change in location-attributed annualised system
cost against the baseline, in bn €/yr. Cost is booked to the location of
the asset that incurs it. The total splits in two. The CDR-chain column
sums DAC and BECCS capture capital and transport and storage, the part
the credit pays for directly. The wider-system column sums renewable
generation, biomass feedstock, and fossil and other effects, none of
which a per-tonne capture credit reaches. Marginal terms on the DAC and
BECCS capture links are excluded throughout, as they carry the
build-time ETS cancellation transfer rather than a resource cost
(Appendix [7.1.3](#subsubsec:app_carbon_accounting){reference-type="ref"
reference="subsubsec:app_carbon_accounting"}). The Captured column is
gross routed to storage from all point sources in the country, in
Mt /yr. It includes fossil and industrial capture, not only credited
CDR. Where industrial capture is present the two diverge: Spain in 2050
captures 50 Mt in total but only 40 Mt as DACCS and BECCS. The Stored
column is geological sequestration in the country, also in Mt /yr.

Cost on the non-spatial European commodity bus, which prices oil and
other primary commodities at a single node, is not attributed to any
country. That bus carries 2.2 bn €/yr in 2050, almost all of it fossil
oil and wider-system, the fourth-largest single entry that year. It is
small in 2030 and 2040 (about 0.1 bn €/yr or less). The listed country
rows therefore do not sum to the European total of about 46 bn €/yr in
2050. The listed rows sum to 44.7 bn €/yr. Adding the commodity bus and
subtracting the omitted net-decrease countries (about 1.0 bn €/yr in
total) recovers the European total.

Countries below the 0.05 bn €/yr reporting threshold are omitted from
the tables. These omitted entries are small by construction and can be
either positive or negative. In 2040 and 2050, several omitted countries
see small net cost decreases relative to the zero-credit baseline. The
credit reorders the least-cost solution: capture, storage, and the
generation feeding them concentrate where they are cheapest, so some
countries build less or import more. The European total still rises,
while the listed rows show where the material positive cost increases
occur.

The same mechanism produces the negative entries within listed
countries. A negative column marks a component whose cost falls below
baseline even where the country total rises. Italy in 2050 is the
clearest case. Its CDR-chain cost rises 0.9 bn €/yr while its
wider-system cost falls 0.5 bn €/yr, as cheaper power displaces local
generation.

Three roles emerge. Hosts such as the United Kingdom and Denmark carry
the largest increases and host most capture and storage, almost all of
it CDR chain, so the credit reaches it. Non-hosting payers such as
Ukraine and Finland add generation for the wider system with no capture
of their own, so the credit never reaches their cost. France sits at the
edge of this group, adding about 0.9 bn €/yr of generation and biomass
cost while hosting about 6 Mt of capture, around 3 Mt of it CDR.
Non-hosting beneficiaries such as Austria reach the target while their
own cost falls, as do most other omitted countries.

::: {#tab:app_country_cost_2030}
+------------+--------------------------------------+-----------------------+
|            | Cost change \[bn €/yr\]              | Hosted \[Mt /yr\]     |
+:===========+==========:+==========:+=============:+==========:+==========:+
| 2-4(lr)5-6 | Total     | CDR chain | Wider system | Captured  | Stored    |
| Country    |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| United     | 0.55      | 0.54      | 0.02         | 29.7      | 29.7      |
| Kingdom    |           |           |              |           |           |
| (GB)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Sweden     | 0.36      | 0.36      | -0.01        | 5.0       | 5.0       |
| (SE)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Italy (IT) | 0.23      | 0.20      | 0.03         | 2.7       | 2.7       |
+------------+-----------+-----------+--------------+-----------+-----------+
| Spain (ES) | 0.22      | 0.22      | 0.01         | 2.9       | 2.9       |
+------------+-----------+-----------+--------------+-----------+-----------+
| Denmark    | 0.15      | 0.13      | 0.03         | 4.3       | 4.3       |
| (DK)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Norway     | 0.11      | 0.12      | 0.00         | 13.0      | 13.0      |
| (NO)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Germany    | 0.09      | 0.08      | 0.02         | 2.3       | 2.3       |
| (DE)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Portugal   | 0.08      | 0.08      | 0.00         | 1.0       | 1.0       |
| (PT)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+

: Country-level system-cost change and hosted carbon removal, medium
cost, 2030. Columns defined in the text.
:::

::: {#tab:app_country_cost_2040}
+------------+--------------------------------------+-----------------------+
|            | Cost change \[bn €/yr\]              | Hosted \[Mt /yr\]     |
+:===========+==========:+==========:+=============:+==========:+==========:+
| 2-4(lr)5-6 | Total     | CDR chain | Wider system | Captured  | Stored    |
| Country    |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Denmark    | 10.51     | 9.91      | 0.59         | 151.7     | 159.6     |
| (DK)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| United     | 5.08      | 3.42      | 1.65         | 67.0      | 66.4      |
| Kingdom    |           |           |              |           |           |
| (GB)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Germany    | 1.71      | 0.59      | 1.12         | 16.1      | 16.6      |
| (DE)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Ireland    | 1.66      | 0.77      | 0.89         | 14.2      | 14.1      |
| (IE)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Spain (ES) | 1.41      | 0.81      | 0.60         | 19.7      | 5.9       |
+------------+-----------+-----------+--------------+-----------+-----------+
| France     | 1.32      | -0.04     | 1.37         | 0.9       | 0.0       |
| (FR)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Sweden     | 0.81      | 0.87      | -0.06        | 16.5      | 8.7       |
| (SE)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Ukraine    | 0.77      | 0.00      | 0.77         | 0.0       | 0.0       |
| (UA)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Poland     | 0.63      | 0.13      | 0.50         | 3.0       | 3.0       |
| (PL)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Finland    | 0.41      | 0.00      | 0.42         | 0.0       | 0.0       |
| (FI)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Greece     | 0.31      | 0.08      | 0.23         | 0.9       | 0.9       |
| (GR)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Austria    | 0.16      | 0.00      | 0.16         | 0.0       | 0.0       |
| (AT)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Belgium    | 0.11      | 0.15      | -0.05        | 2.4       | 0.0       |
| (BE)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+
| Latvia     | 0.08      | 0.11      | -0.03        | 1.8       | 1.8       |
| (LV)       |           |           |              |           |           |
+------------+-----------+-----------+--------------+-----------+-----------+

: Country-level system-cost change and hosted carbon removal, medium
cost, 2040. Columns defined in the text.
:::

::: {#tab:app_country_cost_2050}
+-------------+--------------------------------------+-----------------------+
|             | Cost change \[bn €/yr\]              | Hosted \[Mt /yr\]     |
+:============+==========:+==========:+=============:+==========:+==========:+
| 2-4(lr)5-6  | Total     | CDR chain | Wider system | Captured  | Stored    |
| Country     |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+
| United      | 21.94     | 17.02     | 4.92         | 288.9     | 282.5     |
| Kingdom     |           |           |              |           |           |
| (GB)        |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Denmark     | 12.43     | 12.33     | 0.10         | 197.5     | 197.8     |
| (DK)        |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Spain (ES)  | 4.19      | 1.14      | 3.05         | 50.4      | 8.2       |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Norway (NO) | 1.55      | 1.25      | 0.30         | 27.5      | 27.2      |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Ireland     | 1.39      | 0.78      | 0.61         | 18.7      | 10.7      |
| (IE)        |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+
| France (FR) | 0.92      | 0.03      | 0.89         | 6.4       | 0.0       |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Germany     | 0.55      | 0.56      | -0.01        | 18.3      | 16.3      |
| (DE)        |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Italy (IT)  | 0.47      | 0.93      | -0.46        | 36.1      | 15.0      |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Netherlands | 0.40      | 0.34      | 0.06         | 9.9       | 18.3      |
| (NL)        |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Poland (PL) | 0.23      | 0.04      | 0.19         | 4.6       | 3.0       |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Greece (GR) | 0.23      | 0.14      | 0.09         | 7.1       | 0.9       |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Ukraine     | 0.17      | 0.00      | 0.16         | 0.0       | 0.0       |
| (UA)        |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Sweden (SE) | 0.16      | 0.67      | -0.51        | 10.8      | 10.3      |
+-------------+-----------+-----------+--------------+-----------+-----------+
| Finland     | 0.08      | 0.03      | 0.05         | 2.0       | 0.0       |
| (FI)        |           |           |              |           |           |
+-------------+-----------+-----------+--------------+-----------+-----------+

: Country-level system-cost change and hosted carbon removal, medium
cost, 2050. Columns defined in the text.
:::

### Robustness check of target-binding credit price to the ETS path {#app:sensitivity_ets}

The standalone credit cancels the direct ETS value of a removed tonne
(Section [3.1.1](#subsubsec:supply-modules){reference-type="ref"
reference="subsubsec:supply-modules"}), so the target-binding credit
price that clears the EU CDR target should be largely independent of the
ETS level. This is tested by re-running the medium-cost supply sweep
under the low (70/130/400) and high (160/400/630 €/t) trajectories
alongside the medium case
(Appendix [7.1.12](#subsubsec:app_ets_prices){reference-type="ref"
reference="subsubsec:app_ets_prices"}), reading the target-binding
credit price on the same 50 €/t grid used throughout.

The three trajectories bracket the carbon-price uncertainty
(Figure [16](#fig:app_ets_trajectories){reference-type="ref"
reference="fig:app_ets_trajectories"}). The low path follows the POLES
core projection of Enerdata [@enerdata_carbon_2025], rising from about
75 €/t in 2025 to roughly 130 €/t in 2040, with 2050 set to the 400 €/t
lower bound of that study's carbon-removal sensitivity. The high path is
an upper bracket assembled from several sources. The 2030 value
(160 €/t) takes the high-ambition view of @lseg_eu_2023 under the EU 90%
2040 target, supported by @abn_amro_group_economics_esg_2025 (145 €/t)
and @bnef_eua_2025 (about 149 €/t). The 2040 value (400 €/t) reflects
the same @lseg_eu_2023 statement, and 2050 (630 €/t) is the upper bound
of the Enerdata removals sensitivity. The low and high 2050 values are
the two ends of one Enerdata band.

<figure id="fig:app_ets_trajectories" data-latex-placement="H">
<img
src="./Figures/Appendix/fig_ets_sensitivity_trajectories_style_b.png"
style="width:60.0%" />
<figcaption>EU ETS price trajectories tested in the robustness check:
low (70/130/400), medium (119/279/463), and high (160/400/630) €/t, with
the shaded sensitivity range. The medium path is the main-analysis
trajectory.</figcaption>
</figure>

Two limits of the high path are noted. No modelled trajectory in the
surveyed literature puts the 2050 price above the medium 463 €/t except
the Enerdata upper bound, so the high-path 2050 value is a stress bound,
not a consensus estimate. The same group's scarcity run reaches only
about 400 €/t after 2050 [@pahle_emerging_2025]. The 2040 value of
400 €/t is an ambition statement, against modelled 2040 values near 130
(Enerdata), 180 (CAKE), and 279 €/t (medium). These soft values do not
affect the results, since the target-binding price is 100 €/t on the
high path in every year.

The target-binding credit price is nearly invariant to the ETS path. It
is 100 €/t on every trajectory and horizon except the central path in
2040, which clears one grid step higher at 150 €/t
(Table [33](#tab:app_ets_binding){reference-type="ref"
reference="tab:app_ets_binding"}). The ETS price spans 400 to 630 €/t in
2050 while the binding credit price stays at 100 €/t on all three paths,
far below the carbon price from 2040 onward. Cost-optimal deployment is
set by the credit price, not the carbon price, as intended by the
standalone-credit design.

::: {#tab:app_ets_binding}
  ETS trajectory (2030/2040/2050)    2030   2040   2050
  --------------------------------- ------ ------ ------
  low (70/130/400)                   100    100    100
  medium (119/279/463)               100    150    100
  high (160/400/630)                 100    100    100

  : Target-binding CDR credit price (€/t ) under low, medium, and high
  ETS trajectories, medium-cost sensitivity, on the 50 €/t price grid.
  The target-binding credit price is the lowest credit price at which
  credited CDR meets the EU rCDR targets (22/160/514 Mt /yr).
:::

The ETS trajectory still shapes the wider system, but indirectly. A
higher carbon price raises fossil and process point-source CCS and
raises the biomass opportunity cost that governs the BECCS-to-DACCS
handover (Section [5.1](#sec:disc-trigger){reference-type="ref"
reference="sec:disc-trigger"}). These channels move the deployment and
capture mix, but they leave the credit price needed to reach the removal
goal unchanged, because credited removal takes first claim on geological
storage under the waterfall attribution
(Appendix [7.1.4](#subsubsec:app_waterfall){reference-type="ref"
reference="subsubsec:app_waterfall"}).

The lone exception to this robustness is the non-monotonic 2040
response, where the central path clears one grid step above both the low
and high paths. This reversal lies within the 50 €/t resolution of the
sweep and does not reflect a monotonic price signal. At the 160 Mt goal,
DAC enters as a near-discrete block, so credited deployment crosses the
target within a single price step, and small differences in the wider
system across the three trajectories shift the crossing between the 100
and 150 €/t steps. A finer credit-price grid would resolve the 2040
binding price to a single value and remove the apparent reversal.

Table [34](#tab:app_ets_headline){reference-type="ref"
reference="tab:app_ets_headline"} gives the same headline metrics under
the low, medium, and high ETS paths at medium cost. Each column block
holds one ETS path across the three horizons. The medium block
reproduces the medium-cost case of
Table [22](#tab:app_cost_headline){reference-type="ref"
reference="tab:app_cost_headline"}, the main analysis.

The credit-driven quantities remain nearly flat across the ETS paths.
CDR, physical storage, and the added system cost move little, with the
2050 cost increase holding near 46 to 55 bn €/yr. Fossil & process
capture and the electricity price are the metrics that respond to the
carbon price. Fossil & process capture rises from 8 to 47 Mt /yr in 2030
and from 57 to 86 Mt /yr in 2040, because a higher carbon price makes
point-source capture pay. The average electricity price rises with it,
reaching 110 €/MWh on the high path in 2050 against 54 €/MWh on the low
path. The average heat price moves the same way, reaching 100 €/MWh on
the high path in 2050 against 56 €/MWh on the low path. These are the
indirect channels through which the ETS shapes the system, while the
credit price needed to reach the goal stays almost unchanged.

::: {#tab:app_ets_headline}
+--------------------+--------------------+--------------------+--------------------+
|                    | low                | medium             | high               |
+:===================+:====:+:====:+:====:+:====:+:====:+:====:+:====:+:====:+:====:+
| 2-4(lr)5-7(lr)8-10 | 2030 | 2040 | 2050 | 2030 | 2040 | 2050 | 2030 | 2040 | 2050 |
| Metric             |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Binding credit     | 100  | 100  | 100  | 100  | 150  | 100  | 100  | 100  | 100  |
| price (€/t)        |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Credited CDR       | 45   | 280  | 600  | 23   | 275  | 584  | 23   | 249  | 600  |
| (Mt /yr)           |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| of which DACCS     | 0    | 218  | 539  | 0    | 208  | 551  | 0    | 225  | 576  |
+--------------------+------+------+------+------+------+------+------+------+------+
| of which BECCS     | 45   | 62   | 61   | 23   | 67   | 32   | 23   | 24   | 24   |
+--------------------+------+------+------+------+------+------+------+------+------+
| Fossil & process   | 8    | 57   | 93   | 47   | 61   | 116  | 47   | 86   | 118  |
| CCS (Mt /yr)       |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Average            | 63   | 79   | 54   | 68   | 98   | 69   | 72   | 106  | 110  |
| electricity price  |      |      |      |      |      |      |      |      |      |
| (€/MWh)            |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| Average heat price | 38   | 68   | 56   | 40   | 92   | 52   | 41   | 97   | 100  |
| (€/MWh)            |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+
| $\Delta$ system    | 4    | 22   | 51   | 2    | 23   | 46   | 2    | 23   | 55   |
| cost (bn €/yr)     |      |      |      |      |      |      |      |      |      |
+--------------------+------+------+------+------+------+------+------+------+------+

: Headline metrics at the target-binding credit price by ETS path (low
70/130/400, medium 119/279/463, high 160/400/630 €/t), medium cost.
$\Delta$ system cost excludes the ETS-cancellation transfer, against
each path's own zero-credit case. medium reproduces the medium-cost case
of Table [22](#tab:app_cost_headline){reference-type="ref"
reference="tab:app_cost_headline"}.
:::
