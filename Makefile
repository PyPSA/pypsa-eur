# SPDX-FileCopyrightText: : 2021-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: CC0-1.0

.PHONY: _conda_check install install-fixed test clean-tests reset

# Helper: Check if conda or mamba is installed and set CONDA_OR_MAMBA variable
_conda_check:
	@# Check if conda or mamba is installed and set CONDA_OR_MAMBA variable
	@if command -v conda &> /dev/null; then \
		echo "Conda detected, using Conda..."; \
		$(eval CONDA_OR_MAMBA := conda) \
	elif command -v mamba &> /dev/null; then \
		echo "Conda not found, but Mamba detected. Using Mamba..."; \
		$(eval CONDA_OR_MAMBA := mamba) \
	else \
		echo "Neither Conda nor Mamba is installed. Please install one of them and retry."; \
		exit 1; \
	fi

# Install the environment
install: _conda_check
	@$(CONDA_OR_MAMBA) env create -f envs/environment.yaml
	@$(CONDA_OR_MAMBA) run -n pypsa-eur pre-commit install

# Install fixed environment
install-fixed: _conda_check
	@$(CONDA_OR_MAMBA) env create -f envs/environment.fixed.yaml
	@$(CONDA_OR_MAMBA) run -n pypsa-eur pre-commit install

# Run default tests
test:
	set -e
	snakemake -call solve_elec_networks --configfile config/test/config.electricity.yaml --rerun-triggers=mtime
	snakemake -call all --configfile config/test/config.overnight.yaml --rerun-triggers=mtime
	snakemake -call all --configfile config/test/config.myopic.yaml --rerun-triggers=mtime
	snakemake -call make_summary_perfect --configfile config/test/config.perfect.yaml --rerun-triggers=mtime
	snakemake -call all --configfile config/test/config.scenarios.yaml --rerun-triggers=mtime -n
	echo "All tests completed successfully."

# Cleans all output files from tests
clean-tests:
	snakemake -call solve_elec_networks --configfile config/test/config.electricity.yaml --rerun-triggers=mtime --delete-all-output
	snakemake -call all --configfile config/test/config.overnight.yaml --rerun-triggers=mtime --delete-all-output
	snakemake -call all --configfile config/test/config.myopic.yaml --rerun-triggers=mtime --delete-all-output
	snakemake -call make_summary_perfect --configfile config/test/config.perfect.yaml --rerun-triggers=mtime --delete-all-output
	snakemake -call all --configfile config/test/config.scenarios.yaml --rerun-triggers=mtime -n --delete-all-output

# Removes all created files except for large cutout files (similar to fresh clone)
reset:
	@echo "Do you really wanna continue? This will remove logs, resources, benchmarks, results, and .snakemake directories (y/n): " && \
	read ans && [ $${ans} = y ] && ( \
		rm -r ./logs || true; \
		rm -r ./resources || true; \
		rm -r ./benchmarks || true; \
		rm -r ./results || true; \
		rm -r ./.snakemake || true; \
		echo "Reset completed." \
	) || echo "Reset cancelled."
