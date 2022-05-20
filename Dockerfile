# FROM python: "3.10-alpine"
FROM continuumio/miniconda3
MAINTAINER HectorMSanchezC <sanchez.hmsc@berkeley.edu>
###############################################################################
# Setup Structure
###############################################################################
RUN apt-get update \
    && apt-get install nano \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir MGSurvE \
    && mkdir MGSurvE/Paper \
    && mkdir MGSurvE/Demos
WORKDIR /MGSurvE
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
