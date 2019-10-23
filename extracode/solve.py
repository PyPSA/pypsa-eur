#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 17:07:16 2019

@author: hofmann
"""

import pypsa

solver_options = {'threads': 10, 'method': 2,
                  'crossover': 0,
                  'BarConvTol': 1.e-5,
                  'FeasibilityTol': 1.e-5,
                  'AggFill': 0,
                  'PreDual': 0,
                  'GURO_PAR_BARDENSETHRESH': 200}

for LV in [1.0, 1.09, 1.18, 1.25]:
    n = pypsa.Network(f'/home/vres/data/fabian/pypsa-eur/networks/elec_s_512_lv{LV}_Co2L-3H.nc')
    n.add('GlobalConstraint', 'lv_limit', type='transmission_volume_expansion_limit',
          sense='<=', constant=n.linevolumelimit, carrier_attribute='AC, DC')
