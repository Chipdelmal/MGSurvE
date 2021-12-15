Benchmarks
------------

The size of the landscape along with the number of traps in it determine the total processing power required for the optimization of the traps' placements.

For our tests, we ran a regular square-grid landscape with an increasing number of points and traps. 
We ran this code with the base `DEAP <https://deap.readthedocs.io/en/master/>` implementation of a genetic algorithm with Gaussian mutation, blend crossover and tournament selection.
The algorithm was run for 1000 generations (and five repetitions per scenario), in a server with 88 cores and 3GHz processors (two 22 processor sockets)


The scaling on landscape size follows this behavior:

.. image:: ../../img/bench_PointsVTime.jpg
    :width: 60%
    :align: center

While the scaling on traps number adheres to the following:

.. image:: ../../img/bench_TrapsVTime.jpg
    :width: 60%
    :align: center

The code used for these benchmarks can be found `in this link <https://github.com/Chipdelmal/MGSurvE/tree/main/MGSurvE/benchmarks>`_
