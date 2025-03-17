# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

.ONESHELL:

.PHONY: _conda_check install install-pinned-linux install-pinned-windows install-pinned-macos test checks clean-tests reset

# Helper: Check if conda or mamba is installed and set CONDA_OR_MAMBA variable
# mamba is preferred over conda if both are installed, unless PYPSA_PREFER_CONDA is set to true
_conda_check:
	@# Check prefer_conda environment variable
	@if [ "$$PYPSA_PREFER_CONDA" = "true" ]; then \
		if command -v conda &> /dev/null; then \
			echo "Conda preferred and detected. Using Conda..."; \
			$(eval CONDA_OR_MAMBA := conda) \
		elif command -v mamba &> /dev/null; then \
			echo "Conda preferred but not found. Using Mamba..."; \
			$(eval CONDA_OR_MAMBA := mamba) \
		else \
			echo "Neither Conda nor Mamba is installed. Please install one of them and retry."; \
			exit 1; \
		fi \
	else \
		if command -v mamba &> /dev/null; then \
			echo "Mamba detected, using Mamba..."; \
			$(eval CONDA_OR_MAMBA := mamba) \
		elif command -v conda &> /dev/null; then \
			echo "Mamba not found, but Conda detected. Using Conda..."; \
			$(eval CONDA_OR_MAMBA := conda) \
		else \
			echo "Neither Conda nor Mamba is installed. Please install one of them and retry."; \
			exit 1; \
		fi \
	fi

# Install environment
# E.g. make install or make install name=myenv
install: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/environment.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install
# Install pinned environment
install-pinned-linux: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/linux-pinned.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install
install-pinned-windows: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/windows-pinned.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install
install-pinned-macos: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/macos-pinned.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install


# Run default tests
test:
	set -e
	snakemake solve_elec_networks --configfile config/test/config.electricity.yaml
	snakemake --configfile config/test/config.overnight.yaml
	snakemake --configfile config/test/config.myopic.yaml
	snakemake make_summary_perfect --configfile config/test/config.perfect.yaml
	snakemake --configfile config/test/config.scenarios.yaml -n
	echo "All tests completed successfully."

unit-test:
	pytest test

# Cleans all output files from tests
clean-tests:
	snakemake solve_elec_networks --configfile config/test/config.electricity.yaml --delete-all-output
	snakemake --configfile config/test/config.overnight.yaml --delete-all-output
	snakemake --configfile config/test/config.myopic.yaml --delete-all-output
	snakemake make_summary_perfect --configfile config/test/config.perfect.yaml --delete-all-output
	snakemake --configfile config/test/config.scenarios.yaml -n --delete-all-output

# Removes all created files except for large cutout files (similar to fresh clone)
reset:
	@echo "Do you really wanna continue? This will remove config/config.yaml, logs, resources, benchmarks, results, and .snakemake directories (y/n): " && \
	read ans && [ $${ans} = y ] && ( \
		rm -r ./logs || true; \
		rm -r ./resources || true; \
		rm -r ./benchmarks || true; \
		rm -r ./results || true; \
		rm -r ./.snakemake || true; \
		rm ./config/config.yaml || true; \
		echo "Reset completed." \
	) || echo "Reset cancelled."
