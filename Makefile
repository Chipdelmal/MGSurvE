
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
	- pip install .
	- pip install pytest
	- conda config --add channels conda-forge
	- conda install -c conda-forge deap -y
	- conda install -c conda-forge libpysal -y
	- conda install -c conda-forge cartopy -y

###############################################################################
# Docker
###############################################################################
docker_release:
	- docker build -t chipdelmal/mgsurve:$(version) .
	- docker push chipdelmal/mgsurve:$(version)
	- docker build -t chipdelmal/splatstats:latest .
	- docker push chipdelmal/splatstats:latest

docker_run:
	- docker run \
		-v "$(pwd)"/MGS_sims/paper:/MGSurvE/Paper/sims_out \
		-v "$(pwd)"/MGS_sims/demos:/MGSurvE/Demos/demos_out \
		-it mgsurve:dev bash

docker_build:
	- docker rmi mgsurve:dev -f
	- docker build -t mgsurve:dev .


###############################################################################
# Full Release
###############################################################################
pypi-docker_release:
	- make pypi
	- make dockerExport