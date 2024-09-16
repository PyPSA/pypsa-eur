# SPDX-FileCopyrightText: : 2021-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: CC0-1.0

set -x && \

snakemake -call solve_elec_networks --configfile config/test/config.electricity.yaml --rerun-triggers=mtime && \
snakemake -call all --configfile config/test/config.overnight.yaml --rerun-triggers=mtime && \
snakemake -call all --configfile config/test/config.myopic.yaml --rerun-triggers=mtime && \
snakemake -call make_summary_perfect --configfile config/test/config.perfect.yaml --rerun-triggers=mtime && \
snakemake -call all --configfile config/test/config.scenarios.yaml --rerun-triggers=mtime -n && \

set +x
