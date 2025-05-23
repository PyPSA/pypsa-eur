# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 16:20:06 2025

@author: alice
"""
import pypsa

percorso = r"C:\Users\alice\Desktop\CMCC\pypsa-adb-industry\resources\networks\base_s_39___2030.nc"

n = pypsa.Network(percorso)
n.optimize.create_model()
n.model.solve(solver_name="gurobi")
n.model.compute_infeasibilities()
