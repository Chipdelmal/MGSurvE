Installation
=====

It is strongly recommended to install and use `conda <https://docs.conda.io/en/latest/miniconda.html>`_ for environment management, as some of MGSurvE's dependencies might clash with specific packages versions.
These installation instructions assume conda is already installed in the target system.


Uneventful (tested on Linux)
--------------------------

Create a clean :code:`python=3.9` environment:


.. code-block:: console

   conda create -n MGSurvE python="3.10" -y



Activate the environment, and install MGSurvE using pip:

.. code-block:: console

   conda activate MGSurvE
   pip install MGSurvE


If this ran correctly, try importing the package from the terminal with:

.. code-block:: console

   python
   import MGSurvE as srv


If `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ is not currently installed, we will get a warning that we can safely ignore unless we are planning to work with realistic geo-spatial landscapes. 
In case we are interested in such projects, we can run:


.. code-block:: console

   conda install cartopy -y


If any of these steps resulted in an error, have a look at the next section in this installation guide.


Dependency Errors (sometimes happens on MacOS)
--------------------------


Sometimes, either `DEAP <https://deap.readthedocs.io/en/master/>`_, `libpysal <https://pysal.org/libpysal/>`_, or `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ decide not to play nicely.
For these cases, the following procedure might help.


In a fresh environment on :code:`python=3.9` (or one in which the installation of MGSurvE failed), run the following commands to install `DEAP <https://deap.readthedocs.io/en/master/>`_ and `libpysal <https://pysal.org/libpysal/>`_ (both of which are required for MGSurvE to work):


.. code-block:: console

   conda install deap -y
   conda install libpysal -y


With these packages installed, we can proceed and install MGSurvE:


.. code-block:: console

   pip install MGSurvE

which should be enough to get us started.



Finally, if we are interested in the map-plotting routines, we will need `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ installed. The easiest way to install the dependency is through `anaconda <https://www.anaconda.com/products/individual>`_:


.. code-block:: console

   conda install cartopy -y

   
If this installation fails, or if :code:`import MGSurvE` fails on python due to errors in the dependency, we might need to have a look at their `documentation <https://scitools.org.uk/cartopy/docs/latest/installing.html>`_ for more information.


Most Bulletproof Installation Method
--------------------------

If either of these approaches is failing, try the following chain of commands:


.. code-block:: console

   conda create -n MGSurvE python="3.9"
   conda activate MGSurvE
   conda install deap -y
   pip install MGSurvE


And for the optional dependencies (within the same environment):

.. code-block:: console

   conda install libpysal -y
   conda install cartopy -y