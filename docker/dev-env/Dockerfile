# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

FROM ghcr.io/prefix-dev/pixi:0.79.0

LABEL org.opencontainers.image.source=https://github.com/PyPSA/pypsa-eur

# ca-certificates: pixi >=0.64 uses system CA certs, base image ships none
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates bash git && \
    rm -rf /var/lib/apt/lists/*

# Install pixi env at /opt/pixi-env (allow for bind mounts)
COPY pixi.toml pixi.lock /opt/pixi-env/
RUN cd /opt/pixi-env && pixi install -e default && pixi clean cache --yes

WORKDIR /workspace

# Link pre-installed pixi env into workspace, activate environment and run command
RUN echo '[ ! -e /workspace/.pixi ] && ln -s /opt/pixi-env/.pixi /workspace/.pixi' > /entrypoint.sh && \
    pixi shell-hook -e default --manifest-path /opt/pixi-env/pixi.toml >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]
CMD ["bash"]
