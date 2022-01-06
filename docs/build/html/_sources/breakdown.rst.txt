Pkg Breakdown
=====

`MGSurvE <https://github.com/Chipdelmal/MGSurvE>`_ is built around a `landscape <../html/generated/MGSurvE.landscape.html#module-MGSurvE.landscape>`_
object which stores the information of the sites' and traps' positions, along with the mosquito movement information. 

.. image:: ../../img/MGSurvEDiagSingleSex.jpg


Workflow
------------

A general workflow in `MGSurvE <https://github.com/Chipdelmal/MGSurvE>`_ looks as follows:

1. Define point-set 
    a. Positions (pseudo-random or field data)
    b. Types (pseudo-random or field data)
2. Define mosquito movement type
    a. Movement kernel (package-provided or custom implementation)
3. Define traps'settings
    a. Positions (pseudo-random or field data)
    b. Types (pseudo-random or field data)
    c. Movable/Immovable (for optimization purposes)
4. Instantiate landscape object with info from previous steps
5. Setup `DEAP <https://deap.readthedocs.io/en/master/>`_ for Genetic-Algorithm (GA) optimization
    a. Register individual's and population creator
    b. Register mutation operator
    c. Register crossover operator
    d. Register selector operator
    e. Register fitness evaluator 
6. Run the GA optimizer
7. Update landscape and save results


Components 
------------


The :code:`landscape` object contains several elements that make it easy to define, calculate, plot and export our study sites.

Point and Traps Coordinates (:code:`lnd.pointCoords`, :code:`lnd.trapsCoords`) 
~~~~~~~~~~~~~~~~~~~~~~

These two numpy arrays contain the information on the locations of points of interest in the landscape. Their ordering matters in terms of how secondary objects are calculated and stored throughout the code 
(distances and migration matrices preserve the order defined in these structures).


Kernel Function (:code:`lnd.kernelFunction`)
~~~~~~~~~~~~~~~~~~~~~~

This function defines the probability for individuals to move from one point to another in a time-step. In their most basic form, they calculate a probability :code:`p` as a function of the distance :code:`d`
between points. This :code:`kernel function` is run to calculate the :code:`lnd.migrationMatrix` from the :code:`lnd.distanceMatrix`. It can be defined to take into account other landscape features (such as altitude)
if the parameters are provided (this would not pose a substantial computational burden as it is run only once to setup the landscape).


.. image:: ../../img/expo.png


Traps Kernels (:code:`lnd.trapsKernels`)
~~~~~~~~~~~~~~~~~~~~~~

Additional functions that calculate the relation between distance :code:`d` from the trap to the site from which the individual is moving from, and probability :code:`p` to fall into a trap (depending on the type of trap to account for attractiveness). 
These functions are stored in a dictionary where each entry defines the properties for every trap type. 

More complicated functions can be used to account for complicated features such as wind, elevation, land-type, and more; but, as these functions need to be run on each iteration of optimization cycles, 
the more complicated they become, the more processing power that will be required to compensate for the complexity.


Distances Matrix (:code:`lnd.distanceMatrix`)
~~~~~~~~~~~~~~~~~~~~~~

This numpy array contains the distances between all the points in the landscape in the order that they are stored in the 
:code:`pointCoords`. This matrix is calculated point-wise by using the :code:`distance function` provided to the :code:`landscape` object.


Migration Matrix (:code:`lnd.migrationMatrix`)
~~~~~~~~~~~~~~~~~~~~~~

This matrix contains the probabilities of individuals to migrate from point :code:`a` (row) to point :code:`b` (column) across
the landscape in a time-step. This matrix is internally calculated using the :code:`kernel function` and the distance between sites.


.. image:: ../../img/01.png


Masked Migration Matrix (:code:`lnd.maskedMigration`)
~~~~~~~~~~~~~~~~~~~~~~

Similar to the :code:`migration matrix` but this matrix takes into account the point-types for the probability of movement 
(as provided :code:`traps mask` array). If no :code:`traps mask` is provided, this matrix is equal to the :code:`migration matrix`.


.. image:: ../../img/02.png


Traps Matrix (:code:`lnd.trapsMigration`)
~~~~~~~~~~~~~~~~~~~~~~

Finally, the :code:`traps matrix` contains the probabilities of individuals moving between all the points of the landscape (including 
the traps).

.. image:: ../../img/03.png



Genetic Algorithm 
------------