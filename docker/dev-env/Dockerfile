# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

FROM ghcr.io/prefix-dev/pixi:0.63.2 AS build

LABEL org.opencontainers.image.source=https://github.com/PyPSA/pypsa-eur

RUN pixi self-update

RUN apt-get update && \
    apt-get install -y --no-install-recommends bash git && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /pypsa-eur

COPY pixi.toml pixi.lock ./

RUN pixi install -e default && pixi clean cache --yes
RUN pixi shell-hook -e default > /shell-hook.sh

# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "$@"' >> /shell-hook.sh

# set the entrypoint to the shell-hook script (activate the environment and run the command)
# no more pixi needed in the prod container
ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
CMD ["bash"]
