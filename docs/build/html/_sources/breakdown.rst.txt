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


Distances Matrix
~~~~~~~~~~~~~~~~~~~~~~


Migration Kernel
~~~~~~~~~~~~~~~~~~~~~~

Masked Migration Kernel
~~~~~~~~~~~~~~~~~~~~~~

Traps Matrix
~~~~~~~~~~~~~~~~~~~~~~