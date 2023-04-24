# FROM python: "3.10-alpine"
FROM continuumio/miniconda3
LABEL maintainer="Hector M. Sanchez C. <sanchez.hmsc@berkeley.edu>"
###############################################################################
# Setup Structure
###############################################################################
RUN apt-get update \
    && apt-get install nano \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir MGSurvE
WORKDIR /MGSurvE
###############################################################################
# Copy Requirements and License
###############################################################################
COPY ./conda/requirements.yml . 
COPY LICENSE .
###############################################################################
# Copy Paper and Demo Experiments Files
###############################################################################
COPY ./MGSurvE/ .
###############################################################################
# Install Packages
###############################################################################
RUN conda update -n base -c defaults conda \
    && conda env update --file requirements.yml -n base --prune \
    && rm requirements.yml \
    && conda install deap -y \ 
    && conda install cartopy -y \
    && conda install libpysal -y \
    && pip install MGSurvE

