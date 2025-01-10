# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

FROM condaforge/mambaforge

LABEL org.opencontainers.image.source https://github.com/PyPSA/pypsa-eur

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

RUN apt-get update && apt-get install -y bash git make

RUN conda --version

WORKDIR /pypsa-eur

COPY ./envs ./temp

RUN conda env create -n pypsa-eur -f temp/linux-pinned.yaml
RUN conda init bash
RUN echo "conda activate pypsa-eur" >> ~/.bashrc

SHELL ["/bin/bash", "--login", "-c"]
ENV PATH=/opt/conda/envs/pypsa-eur/bin:$PATH

RUN rm -r temp
RUN conda clean -afy && \
    rm -rf /tmp/*

CMD ["bash"]
