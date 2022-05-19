# FROM python: "3.10-alpine"
FROM continuumio/miniconda3
###############################################################################
# Setup Structure
###############################################################################
RUN apt-get update \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir MGSurvE \
    && mkdir MGSurvE/Paper \
    && mkdir MGSurvE/Paper/sims_out \
    && mkdir MGSurvE/Demos \
    && mkdir MGSurvE/Demos/demos_out \
    && mkdir MGSurvE/sims_out
WORKDIR /MGSurvE
RUN ln -s ./Paper/sims_out ./sims_out/paper \
    && ln -s ./Demos/demos_out ./sims_out/demos
###############################################################################
# Copy Requirements and License
###############################################################################
COPY ./conda/requirements.yml . 
COPY LICENSE .
###############################################################################
# Copy Paper and Demo Experiments Files
###############################################################################
COPY ./MGSurvE/demos/Paper/GEO ./Paper/
COPY ./MGSurvE/demos/Paper/*.py ./Paper/
COPY ./MGSurvE/demos/Paper/*.sh ./Paper/
COPY ./MGSurvE/demos/*.sh ./Demos/
COPY ./MGSurvE/demos/*.py ./Demos/
###############################################################################
# Install Packages
###############################################################################
RUN conda env update --file requirements.yml -n base --prune \
    && rm requirements.yml
###############################################################################
# Scratch
###############################################################################
# docker run -v "$(pwd)":/MGSurvE/Paper/sims_out -it mgsurve:latest bash