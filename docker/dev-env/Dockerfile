# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

FROM ghcr.io/prefix-dev/pixi:0.84.2

LABEL org.opencontainers.image.source=https://github.com/PyPSA/pypsa-eur

# ca-certificates: pixi >=0.64 uses system CA certs, base image ships none
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates bash git && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Install the pixi env at its runtime path, since pixi envs are not relocatable,
# then move it to /opt/pixi-env so that a bind-mounted /workspace can link it back
COPY pixi.toml pixi.lock ./
RUN pixi install -e default && pixi clean cache --yes && \
    pixi shell-hook -e default > /tmp/shell-hook.sh && \
    mkdir /opt/pixi-env && mv .pixi /opt/pixi-env/

# Link pre-installed pixi env into workspace, activate environment and run command
RUN echo '[ ! -e /workspace/.pixi ] && ln -s /opt/pixi-env/.pixi /workspace/.pixi' > /entrypoint.sh && \
    cat /tmp/shell-hook.sh >> /entrypoint.sh && rm /tmp/shell-hook.sh && \
    echo 'exec "$@"' >> /entrypoint.sh

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]
CMD ["bash"]
