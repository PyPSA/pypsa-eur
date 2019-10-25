# Builder Stage
FROM debian:stretch as builder

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/conda/bin:$PATH

RUN apt-get update --fix-missing && apt-get install -yq --no-install-recommends curl bzip2 xz-utils git ca-certificates

# download pypsa-eur and data dependencies

# while testing use travis feature branch
# RUN git clone https://github.com/PyPSA/pypsa-eur.git /home
RUN git clone --single-branch --branch travis https://github.com/PyPSA/pypsa-eur.git /home

RUN curl -L "https://zenodo.org/record/3517921/files/pypsa-eur-tutorial-data-bundle.tar.xz" -o "/tmp/pypsa-eur-data-bundle.tar.xz"

RUN curl -L "https://zenodo.org/record/3518020/files/pypsa-eur-tutorial-cutouts.tar.xz" -o "/tmp/pypsa-eur-cutouts.tar.xz"

RUN curl -L "https://zenodo.org/record/3518215/files/natura.tiff" -o "/home/resources/natura.tiff"

RUN tar xJf /tmp/pypsa-eur-data-bundle.tar.xz -C /home/data 

RUN tar xJf /tmp/pypsa-eur-cutouts.tar.xz -C /home

# install and set-up conda

RUN curl -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh \
  && bash miniconda.sh -b -p /opt/conda \
  && . /opt/conda/etc/profile.d/conda.sh \
  && conda install python=3.6 conda-build \
  && conda config --add create_default_packages python=3.6 \
  && conda config --add create_default_packages ipython \
  && conda config --add create_default_packages blas=*=openblas \
  && conda env create -f /home/environment.yaml \
  && conda develop $GUROBI_HOME/lib/python3.6_utf32 \
  && conda clean --all -y \
  && find /opt/conda \( -type d -a -name test -o -name tests \) \
      -o \( -type f -a -name '*.pyc' -o -name '*.pyo' \) \
      -exec rm -rf '{}' + \
  && rm -rf /opt/conda/pkgs/*

# install open-source solvers

RUN apt-get install coinor-cbc && conda install -c conda-forge ipopt 

# Final Stage
FROM debian:stretch
# WARNING: Dont use ARG in the final Stage!
# All ARG Values in this stage will be saved in the final image.

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 \
    PATH=/opt/conda/bin:$PATH

COPY --from=builder /opt /opt
COPY --from=builder /home /home

