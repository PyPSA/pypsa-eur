FROM condaforge/mambaforge

RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

WORKDIR /pypsa-eur

COPY ./envs ./envs

RUN conda env create --file envs/environment.yaml

RUN rm -r envs

RUN echo "source activate pypsa-eur" > ~/.bashrc
ENV PATH /opt/conda/envs/pypsa-eur/bin:$PATH

ENTRYPOINT ["tail", "-f", "/dev/null"]
