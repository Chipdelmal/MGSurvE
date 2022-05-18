# FROM python:"3.9-alpine"
FROM continuumio/miniconda3

# RUN addgroup --system mgsurve
# RUN adduser --system user --ingroup mgsurve

WORKDIR /MGSurvE
COPY ./requirements.yml .

# USER root
RUN apt-get update \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir MGSurvE \
    && conda env update --file requirements.yml -n base --prune
# USER user


# ENTRYPOINT ["conda", "run", "-n", "MGSurvE"]

# COPY ./MGSurvE .
# CMD [ "conda", "activate", "MGSurvE" ] 

# RUN pip install . --use-feature=in-tree-build

# RUN conda env update --file ./conda/requirements.yml --name base
# RUN conda config --add channels conda-forge && \
#     conda install -c conda-forge deap -y && \
#     conda install -c conda-forge libpysal -y && \
#     conda install -c conda-forge cartopy -y
# RUN pip install MGSurvE

# docker image ls
# docker run -it mgsurve:latest bash
# docker rm $(docker ps -a -q)
# docker rmi -f $(docker images -q)