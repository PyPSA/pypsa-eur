FROM condaforge/mambaforge

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

WORKDIR /pypsa-eur

COPY . .

RUN conda env create --file envs/environment.yaml

RUN echo "source activate pypsa-eur" > ~/.bashrc
ENV PATH /opt/conda/envs/pypsa-eur/bin:$PATH

# ENTRYPOINT [ "snakemake","--cores","1","solve_networks" ]
# ENTRYPOINT ["tail", "-f", "/dev/null"]
ENTRYPOINT ["tail", "-f", "/dev/null"]
