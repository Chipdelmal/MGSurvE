Installation
=====


To use MGSurvE, install it using pip:

.. code-block:: console

   $ pip install MGSurvE

MGSurvE needs `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ installed. The easiest way to install the dependency is through `anaconda <https://www.anaconda.com/products/individual>`_:


.. code-block:: console

   $ conda install -c conda-forge cartopy

   

Alternatively, dependencies can be installed through our conda environment through its YML:

.. code-block:: console

   $ conda env create -f REQUIREMENTS.yml


or TXT definitions:

.. code-block:: console

   $ conda create -n new MGSurvE --file REQUIREMENTS.txt
