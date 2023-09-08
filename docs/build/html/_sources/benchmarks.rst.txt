Benchmarks
------------

MGSurvE's optimization runtime depends on several factors but the main ones are the number of sites and traps in our landscape.
This is because they are both related to the size of the migration matrix and the number of operations that need to be performed on it. 
Specifically, these variables both impact the calculation of Markov's Fundamental Matrix, which is the slowest part of the optimization process.

To provide a reference on how much time we could expect optimization we provide a `benchmarking routine <https://github.com/Chipdelmal/MGSurvE/blob/main/MGSurvE/benchmarks/dst_benchmark.py>`_.
We ran a set of experiments limiting our server to use 16 cores and ran 5 repetitions of a set of sites-traps combinations (Latin Hypercube Sampling schemes were used) to test how the optimizations time scales (GA parameters were set to :code:`auto` for consistency).
Both optimization routines (discrete and continuous) were ran over 500 generations, as that gets us a good sample of the timing for each parameter combination (after the optimization algorithm's memory allocation and initial function calls have taken place).
We are, however, scaling the timing to 1000 generations (doubling it) because that's closer to the expected number of generations we would run for an optimization task.


.. image:: ../../img/timings_DSC.jpg
    :width: 100%


.. .. |pic2| image:: ../../img/timings_CNT.png
..     :width: 100%


Our server has 44 physical cores (two "Intel(R) Xeon(R) CPU E5-2696 v4 @ 2.20GHz") and 256Gb RAM but, as was previously mentioned, was limited to 16 cores for benchmarking purposes.



.. Distance Functions (Author: Elijah Bartolome)
.. ~~~~~~~~~~~~~~~~~~~~~~

.. 1,000,000 pairs of random points were created. Each point had a valid longitude and latitude value (the longitude was a random number between and while the latitude was a random number between and)

.. The distance between each pair of random points was calculated with each distance function. Each function used the same pairs of points. Here are the times it took to calculate the distance between all 1,000,000 pairs:

.. * Vincenty: 13.3743638 seconds
.. * `Cheap Ruler <https://github.com/mapbox/cheap-ruler>`_: 1.6893626000000026 seconds
.. * Haversine: 2.408093000000001 seconds

.. Here are violin plots of the distribution of execution times for all 1,000,000 pairs of points for each distance function: 

.. .. image:: ../../img/distancePlots5.jpg
..     :width: 100%
..     :align: center


.. Both Haversine and `Cheap Ruler <https://github.com/mapbox/cheap-ruler>`_ are about 10 seconds faster than Vincenty when calculating the 1,000,000 points. So for optimizing runtime, Haversine and `Cheap Ruler <https://github.com/mapbox/cheap-ruler>`_  are ideal with `Cheap Ruler <https://github.com/mapbox/cheap-ruler>`_  having a slight advantage.
.. `Cheap Ruler <https://github.com/mapbox/cheap-ruler>`_ , however, has a wide, problematic range of error. Haversine, in comparison, has an error range orders of magnitude smaller than `Cheap Ruler <https://github.com/mapbox/cheap-ruler>`_ .
.. If one wants to try to minimize runtimes while trying to preserve the accuracy of the distance function, then Haversine is the optimal distance function compared to Vincenty and Haversine.


.. .. image:: ../../img/errorPlots3.jpg
..     :width: 100%
..     :align: center