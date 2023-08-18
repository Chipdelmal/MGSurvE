
SHELL=bash
python=python
pip=pip
tests=.
version:=$(shell $(python) version.py)
sdist_name:=MGSurvE-$(version).tar.gz


###############################################################################
# Unit Tests
###############################################################################
test:
	@$(python) -m pytest -vv $(tests) --disable-pytest-warnings

check_build_reqs:
	@$(python) -c 'import pytest' \
                || ( printf "$(redpip)Build requirements are missing. Run 'make prepare' to install them.$(normal)" ; false )

###############################################################################
# Pypi
###############################################################################
pypi: clean clean_sdist
	set -x \
	&& $(python) setup.py sdist bdist_wheel \
	&& twine check dist/* \
	&& twine upload dist/*

clean_pypi:
	- rm -rf build/

conda_export:
	- pip freeze > ./conda/requirements.txt
	- conda env export | cut -f 1 -d '=' | grep -v "prefix" > ./conda/requirements.yml

conda_update:
	- conda update --all -y
	- pip freeze > ./conda/requirements.txt
	- conda env export | cut -f 1 -d '=' | grep -v "prefix" > ./conda/requirements.yml

###############################################################################
# Docs
###############################################################################
doc:
	- pip install .
	- sphinx-apidoc -f -o docs/source MGSurvE
	- sphinx-build -b html docs/source/ docs/build/html

###############################################################################
# Dev
###############################################################################
develop:
	$(pip) install -e .

clean_develop:
	- $(pip) uninstall -y MGSurvE
	- rm -rf *.egg-info

clean_sdist:
	- rm -rf dist

clean: 
	- clean_develop clean_pypi

dev: 
	- make clean
	- make develop
	- make test

dev_full: 
	- yes | pip install pytest sphinx sphinx_rtd_theme chardet twine wheel jupyter
	- conda update -n base -c defaults conda -y
	- conda update -n base -c conda-forge conda -y
	- conda install -n base conda-libmamba-solver -y
	- conda config --set solver libmamba -y
	- conda config --add channels bioconda -y
	- conda config --add channels conda-forge -y
	- conda install gdal fiona pyproj cartopy libpysal -y
	- conda install -c conda-forge deap nodejs osmnx basemap-data-hires -y
	- pip install .


###############################################################################
# Docker
###############################################################################
docker_release:
	- docker buildx build . \
		--platform=linux/amd64,linux/arm64 \
		-t chipdelmal/mgsurve:$(version) \
		-t chipdelmal/mgsurve:latest \
		--push

docker_run:
	- docker run -it mgsurve:dev bash

docker_build:
	- docker rmi mgsurve:dev -f
	- docker build -t mgsurve:dev .

###############################################################################
# Full Release
###############################################################################
pypi-docker_release:
	- make pypi
	- make docker_release