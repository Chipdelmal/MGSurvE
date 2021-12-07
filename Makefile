define help

Supported targets: prepare, develop, sdist, clean, test, and pypi.

The 'prepare' target installs this project's build requirements into the current virtualenv.

The 'develop' target creates an editable install of this project and its runtime requirements in the
current virtualenv. The install is called 'editable' because changes to the source code
immediately affect the virtualenv.

The 'clean' target undoes the effect of 'develop'.

The 'pypi' target publishes the current commit of this project to PyPI after enforcing that the working
copy and the index are clean, and tagging it as an unstable .dev build.

endef
export help
help:
	@printf "$$help"

SHELL=bash
python=python
pip=pip
tests=.
version:=$(shell $(python) version.py)
sdist_name:=MGSurvE-$(version).tar.gz

develop:
	$(pip) install -e .

clean_develop:
	- $(pip) uninstall -y MGSurvE
	- rm -rf *.egg-info

clean_sdist:
	- rm -rf dist

clean: clean_develop clean_pypi

test:
	@$(python) -m pytest -vv $(tests)

check_build_reqs:
	@$(python) -c 'import pytest' \
                || ( printf "$(redpip)Build requirements are missing. Run 'make prepare' to install them.$(normal)" ; false )

pypi: clean clean_sdist
	set -x \
	&& $(python) setup.py sdist bdist_wheel \
	&& twine check dist/* \
	&& twine upload dist/*

clean_pypi:
	- rm -rf build/

doc:
	- pip install .
	- sphinx-apidoc -f -o docs/source MGSurvE
	- sphinx-build -b html docs/source/ docs/build/html

dev: 
	- make clean
	- make develop
	- make test
	- pip install .

conda:
	- conda list -e > /conda/REQUIREMENTS.txt
	- conda env export > /conda/REQUIREMENTS.yml