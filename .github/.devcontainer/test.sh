#!/bin/bash

ls

pip install git+https://github.com/open-energy-transition/linopy.git@only-generate-problem-files --no-deps

conda install -c conda-forge time

cp -r benchmarks/pypsa/* pypsa-eur/config

mv pypsa-eur/config/solver_benchmark_pypsa_eur.py pypsa-eur/

sed -i '/cf_solving = solving\["options"\]/a\ kwargs["keep_files"] = cf_solving.get("keep_files", True)' pypsa-eur/scripts/solve_network.py

cd pypsa-eur

python solver_benchmark_pypsa_eur.py --configfile config/pypsa-eur-elec-10-lvopt-3h.yaml+