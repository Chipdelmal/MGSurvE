Installation
=====


To use MGSurvE, install it using pip:

.. code-block:: console

   $ pip install MGSurvE

MGSurvE needs `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ installed. The easiest way to install the dependency is through `anaconda <https://www.anaconda.com/products/individual>`_:


.. code-block:: console

   $ conda install cartopy -y

   

Installing `MGSurvE <https://github.com/Chipdelmal/MGSurvE>`_ in a `conda <https://docs.conda.io/en/latest/miniconda.html>`_ environment is strongly suggested, so once it's downloaded, simply run:

.. code-block:: python

    conda create -n MGSurvE python="3.10" -y
    conda activate MGSurvE
    conda install deap -y
    conda install cartopy -y
    pip install MGSurvE

Alternatively, dependencies can be installed through our conda environment through its YML:

.. code-block:: console

   $ conda env create -f REQUIREMENTS.yml


or TXT definitions:

.. code-block:: console

   $ conda create -n new MGSurvE --file REQUIREMENTS.txt



    