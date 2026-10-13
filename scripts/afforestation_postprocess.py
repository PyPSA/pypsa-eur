# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Post-processing redistribution for afforestation carbon dioxide removal links/stores.

Afforestation is solved with free dispatch (fixed p_nom = peak rate,
p_min_pu=0): the LP solver freely distributes the annual CO₂ removal
across the year without a seasonal signal. The solved dispatch is
physically valid (annual total correct, instantaneous rate bounded) but
its temporal shape is arbitrary (LP degeneracy).

This module redistributes stores_t.e and links_t.p0/p1 proportionally
to the precomputed seasonal weights so that temporal plots are
physically meaningful, without re-solving.

Usage
-----
Call once right after loading the solved network, before any temporal
plotting::

    import pypsa
    from scripts.afforestation_postprocess import redistribute_afforestation_seasonal

    n = pypsa.Network("results/.../network_s_90_ec_l_2050.nc")
    seasonal_csv = "resources/afforestation_seasonal_profile_s_90.csv"
    redistribute_afforestation_seasonal(n, seasonal_csv)  # modifies in-place

Only components with carrier == "co2 afforestation" are touched.
All energy-system quantities are unchanged.

Mathematical definition
-----------------------
Given:
  w[t, node]  — seasonal weight from profile CSV (raw, not normalised)
  E[node]     — optimised annual CO₂ store size [tCO₂]  (= e_nom_opt)
  η           — link efficiency (1: the CRCF factor is applied to the potential)

  Δe[t, node] = w[t, node] / Σ_t w[t, node]  ×  E[node]   [tCO₂ per snapshot]
  e[t, node]  = cumsum_t Δe                                  [tCO₂]
  p1[t, node] = −Δe[t, node]                                 [tCO₂/h, bus1 sign]
  p0[t, node] =  Δe[t, node] / η                             [tCO₂/h, bus0 sign]

p0 > 0  → CO₂ withdrawn from "co2 atmosphere"
p1 < 0  → CO₂ injected into "{node} co2 afforestation" bus (PyPSA link convention)
"""

import logging

import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)


def redistribute_afforestation_seasonal(
    n: pypsa.Network,
    seasonal_profile_csv: str,
) -> pypsa.Network:
    """
    Redistribute afforestation CO₂ flows according to seasonal weights.

    Modifies n in-place. Returns n for convenience.

    Parameters
    ----------
    n : pypsa.Network
        Solved network containing afforestation links/stores.
    seasonal_profile_csv : str or Path
        Path to the seasonal profile CSV produced by build_afforestation_potentials.
        Expected columns: node names like "DE0 afforestation".
        Expected index name: "snapshot" (datetime strings).
    """
    affo_links = n.links.index[n.links.carrier == "co2 afforestation"]

    if affo_links.empty:
        logger.warning(
            "No links with carrier 'co2 afforestation' found — "
            "redistribute_afforestation_seasonal is a no-op."
        )
        return n

    # ── Load and align seasonal profile ───────────────────────────────────────
    profile = pd.read_csv(seasonal_profile_csv, index_col="snapshot", parse_dates=True)

    dates = (
        n.snapshots.get_level_values(1)
        if isinstance(n.snapshots, pd.MultiIndex)
        else n.snapshots
    )
    profile_aligned = profile.reindex(dates)

    # ── Ensure time-series frames have columns for all afforestation components ─
    affo_stores = n.stores.index[n.stores.carrier == "co2 afforestation"]
    for col in affo_links:
        if col not in n.links_t.p0.columns:
            n.links_t.p0[col] = 0.0
        if col not in n.links_t.p1.columns:
            n.links_t.p1[col] = 0.0
    for col in affo_stores:
        if col not in n.stores_t.e.columns:
            n.stores_t.e[col] = 0.0

    # ── Redistribute node by node ──────────────────────────────────────────────
    n_done = 0
    for link_name in affo_links:
        node = link_name.removesuffix(" afforestation")
        store_name = f"{node} co2 afforestation"
        profile_col = f"{node} afforestation"

        if store_name not in affo_stores:
            logger.debug(f"No store for link {link_name!r}, skipping.")
            continue

        e_nom_opt = n.stores.at[store_name, "e_nom_opt"]
        if not np.isfinite(e_nom_opt) or e_nom_opt <= 0:
            logger.debug(f"e_nom_opt={e_nom_opt} for {store_name!r}, skipping.")
            continue

        if profile_col not in profile_aligned.columns:
            logger.warning(
                f"Column {profile_col!r} missing from seasonal profile — skipping {node!r}."
            )
            continue

        weights = profile_aligned[profile_col].set_axis(n.snapshots)
        w_sum = weights.sum()
        if w_sum <= 0:
            logger.warning(f"Zero total weight for {node!r}, skipping.")
            continue

        eff = n.links.at[link_name, "efficiency"]
        if not (np.isfinite(eff) and eff > 0):
            eff = 1.0

        # CO₂ stored per snapshot [tCO₂]
        delta_e = weights / w_sum * e_nom_opt

        n.stores_t.e[store_name] = delta_e.cumsum().values
        n.links_t.p0[link_name] = (delta_e / eff).values  # withdrawal from atmosphere
        n.links_t.p1[link_name] = (-delta_e).values  # injection into store bus

        n_done += 1

    logger.info(
        f"redistribute_afforestation_seasonal: redistributed {n_done}/{len(affo_links)} nodes."
    )
    return n
