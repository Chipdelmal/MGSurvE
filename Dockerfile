# FROM python: "3.10-alpine"
FROM continuumio/miniconda3
###############################################################################
# Setup Structure
###############################################################################
RUN apt-get update \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir MGSurvE \
    && mkdir MGSurvE/Paper \
    && mkdir MGSurvE/Paper/sims_out
WORKDIR /MGSurvE
###############################################################################
# Copy Requirements and License
###############################################################################
COPY ./conda/requirements.yml .
COPY LICENSE .
###############################################################################
# Copy Paper Experiments Files
###############################################################################
COPY ./MGSurvE/demos/Paper/*.py ./Paper/
COPY ./MGSurvE/demos/Paper/*.sh ./Paper/
###############################################################################
# Install Packages
###############################################################################
RUN conda env update --file requirements.yml -n base --prune

# RUN pip install . --use-feature=in-tree-build

# RUN conda env update --file ./conda/requirements.yml --name base
# RUN conda config --add channels conda-forge && \
#     conda install -c conda-forge deap -y && \
#     conda install -c conda-forge libpysal -y && \
#     conda install -c conda-forge cartopy -y
# RUN pip install MGSurvE

# docker run -v "$(pwd)":/MGSurvE/Paper/sims_out -it mgsurve:latest bash