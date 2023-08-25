Installation
=====

It is strongly recommended to install and use `conda <https://docs.conda.io/en/latest/miniconda.html>`_ for environment management, as some of MGSurvE's dependencies might clash with specific packages versions.
These installation instructions assume `conda <https://docs.conda.io/en/latest/miniconda.html>`_ is already installed in the target system.


Local Installation
^^^^^^^^^^^^^^^^^^^^^

MGSurvE makes use of some external dependencies to make the most of its capabilities. Some of these, such as `DEAP <https://deap.readthedocs.io/en/master/>`_, are strict requirements but others are optional, such as:

*  `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_: needed for geographical landscapes plots and manipulations (not needed if we will be working with Euclidean landscapes but highly recommended nevertheless)
*  `libpysal <https://pysal.org/libpysal/>`_: not strictly required by MGSurvE but useful if we wanted to do more complex point-sampling operations
*  `osmnx <https://osmnx.readthedocs.io/en/stable/>`_: not required but pretty useful for interacting with `OpenStreetMap <https://www.openstreetmap.org/>`_ data


Taking this into account, we provide some alternatives we have tested to get MGSurvE setup. The core software is installed by following the 
`Uneventful Base Installation <https://chipdelmal.github.io/MGSurvE/build/html/installation.html##uneventful-base-installation>`_, whereas installing the package with all the dependencies is achieved through the 
`Additional Dependencies and Installation Mishaps <https://chipdelmal.github.io/MGSurvE/build/html/installation.html#additional-dependencies-and-installation-mishaps>`_ or 
`Bulletproof Installation Method <https://chipdelmal.github.io/MGSurvE/build/html/installation.html#bulletproof-installation-method>`_, as described in this webpage. Finally, if everything described before fails, we can download a 
`Docker Image <https://hub.docker.com/repository/docker/chipdelmal/mgsurve>`_ following the `docker <https://chipdelmal.github.io/MGSurvE/build/html/installation.html#docker>`_ set of instructions.


Uneventful Base Installation
--------------------------

Create a clean :code:`python=3.9` environment:


.. code-block:: console

   conda create -n MGSurvE python="3.9" -y



Activate the environment, and install MGSurvE using pip:

.. code-block:: console

   conda activate MGSurvE
   pip install MGSurvE


If this ran correctly, try importing the package from the terminal with:

.. code-block:: console

   python
   import MGSurvE as srv


If `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ or `libpysal <https://pysal.org/libpysal/>`_ are not currently installed, we will get a warning that we can safely ignore (see the next section for more info); but if any of these steps resulted in an error, let's have a look at the next section in this installation guide.


Additional Dependencies, and Installation Mishaps
--------------------------


Sometimes, either `DEAP <https://deap.readthedocs.io/en/master/>`_, `libpysal <https://pysal.org/libpysal/>`_, or `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ decide not to play nicely.
For these cases, the following procedure might help.


In a fresh environment on :code:`python=3.9` (or one in which the installation of MGSurvE failed), run the following commands to install `DEAP <https://deap.readthedocs.io/en/master/>`_:


.. code-block:: console

   conda install deap -y


Now, let's install `cartopy <https://scitools.org.uk/cartopy/docs/latest/index.html>`_ (for geo-features plotting), and `libpysal <https://pysal.org/libpysal/>`_ (to generate poisson-distributed pointsets). 
The easiest way to install these dependencies is through `anaconda <https://www.anaconda.com/products/individual>`_:


.. code-block:: console

   conda install cartopy -y
   conda install libpysal -y


With these packages installed, we can proceed and install MGSurvE:


.. code-block:: console

   pip install MGSurvE

which should be enough to get us started. If this installation fails, or if :code:`import MGSurvE` fails on python, we might need to have a look at the next section.


Bulletproof Installation Method
--------------------------

If either of these approaches is failing, try the following chain of commands:


.. code-block:: console

   conda create -n MGSurvE python="3.10"
   conda activate MGSurvE
   conda update -n base -c defaults conda
   conda config --add channels conda-forge 
   conda config --add channels bioconda 
   conda install -n base conda-libmamba-solver 
   conda config --set solver libmamba 
   conda install gdal fiona pyproj cartopy libpysal -y 
   conda install -c conda-forge deap nodejs osmnx basemap-data-hires -y 
   pip install MGSurvE


In case this method still fails, please have a look at the installation instructions on: `DEAP <https://deap.readthedocs.io/en/master/installation.html>`_, `cartopy <https://scitools.org.uk/cartopy/docs/latest/installing.html>`_, and `libpysal <https://pysal.org/libpysal/installation.html>`_; before installing MGSurvE.


Docker
^^^^^^^^^^^^^^^^^^^^^

An alternative way to get `MGSurvE <https://github.com/Chipdelmal/MGSurvE>`_ running is to download the latest `Docker image <https://hub.docker.com/r/chipdelmal/mgsurve>`_.

To get started with this approach, `install Docker <https://docs.docker.com/get-docker/>`_ on your system. Then, go to our `Docker image hub <https://hub.docker.com/r/chipdelmal/mgsurve>`_ and get the latest tag to pull from it. For example:

.. code-block:: console

   docker pull chipdelmal/mgsurve


And then run it in interactive mode:

.. code-block:: console

   docker run -it chipdelmal/mgsurve bash


We can check that `MGSurvE <https://github.com/Chipdelmal/MGSurvE>`_ is correctly installed with its dependencies by running:

.. code-block:: console

   python
   import MGSurvE

Which should run without any error or warnings. To quit our python session, we simply type :code:`exit()` into the interpreter, and we exit the docker promt by typing :code:`exit`.