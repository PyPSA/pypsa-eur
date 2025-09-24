# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

.ONESHELL:

.PHONY: _conda_check install install-pinned-linux install-pinned-windows install-pinned-macos install-lock-linux64 install-lock-linux-arm install-lock-windows install-lock-macos64 install-lock-macos-arm test checks clean-tests reset

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
# Install from conda-lock files (recommended)
install-lock-linux64: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/linux-64.lock.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install

install-lock-windows: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/win-64.lock.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install

install-lock-macos64: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/osx-64.lock.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install

install-lock-macos-arm: _conda_check
	$(CONDA_OR_MAMBA) env create -f envs/osx-arm64.lock.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install

# Install pinned environment (deprecated, but maintained for backward compatibility)
install-pinned-linux: _conda_check
	@echo "WARNING: The 'install-pinned-linux' target is deprecated and will be removed in a future release."
	@echo "Please use 'install-lock-linux64' (for x86_64) or 'install-lock-linux-arm' (for ARM) instead"
	@echo "Or directly use: conda env create -f envs/linux-64.lock.yaml"
	$(CONDA_OR_MAMBA) env create -f envs/linux-64.lock.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install

install-pinned-windows: _conda_check
	@echo "WARNING: The 'install-pinned-windows' target is deprecated and will be removed in a future release."
	@echo "Please use 'install-lock-windows' instead"
	@echo "Or directly use: conda env create -f envs/win-64.lock.yaml"
	$(CONDA_OR_MAMBA) env create -f envs/win-64.lock.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install

install-pinned-macos: _conda_check
	@echo "WARNING: The 'install-pinned-macos' target is deprecated and will be removed in a future release."
	@echo "Please use 'install-lock-macos64' (for Intel) or 'install-lock-macos-arm' (for Apple Silicon) instead"
	@echo "Or directly use: conda env create -f envs/osx-64.lock.yaml"
	$(CONDA_OR_MAMBA) env create -f envs/osx-64.lock.yaml -n $(or $(name), pypsa-eur)
	$(CONDA_OR_MAMBA) run -n $(or $(name), pypsa-eur) pre-commit install


# Run default tests
test:
	set -e
	snakemake -call solve_elec_networks --configfile config/test/config.electricity.yaml
	snakemake -call --configfile config/test/config.overnight.yaml
	snakemake -call --configfile config/test/config.myopic.yaml
	snakemake -call make_summary_perfect --configfile config/test/config.perfect.yaml
	snakemake -call resources/test-elec-clusters/networks/base_s_adm.nc --configfile config/test/config.clusters.yaml
	snakemake -call --configfile config/test/config.scenarios.yaml -n
	snakemake -call plot_power_networks_clustered --configfile config/test/config.tyndp.yaml
	echo "All tests completed successfully."

unit-test:
	pytest test

# Cleans all output files from tests
clean-tests:
	snakemake -call solve_elec_networks --configfile config/test/config.electricity.yaml --delete-all-output
	snakemake -call --configfile config/test/config.overnight.yaml --delete-all-output
	snakemake -call --configfile config/test/config.myopic.yaml --delete-all-output
	snakemake -call make_summary_perfect --configfile config/test/config.perfect.yaml --delete-all-output
	snakemake -call resources/test-elec-clusters/networks/base_s_adm.nc --configfile config/test/config.clusters.yaml --delete-all-output
	snakemake -call --configfile config/test/config.scenarios.yaml -n --delete-all-output
	snakemake -call plot_power_networks_clustered --configfile config/test/config.tyndp.yaml --delete-all-output

# Removes all created files except for large cutout files (similar to fresh clone)
reset:
	@echo "Do you really wanna continue? This will remove logs, resources, benchmarks, results, and .snakemake directories (config/config.yaml will not be deleted) (y/n): " && \
	read ans && [ $${ans} = y ] && ( \
		rm -r ./logs || true; \
		rm -r ./resources || true; \
		rm -r ./benchmarks || true; \
		rm -r ./results || true; \
		rm -r ./.snakemake || true; \
		echo "Reset completed." \
	) || echo "Reset cancelled."
