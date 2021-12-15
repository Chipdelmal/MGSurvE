Landscape Update
------------

Updating the number, position and types of traps is easy. 
Just generate a new dataframe with the traps' information and call the :code:`lnd.updateTraps` function:

.. code-block:: python

    traps = pd.DataFrame({
        'x': [0.5, 3.0, 2.0], 
        'y': [0.0, 0.0, 2.0], 
        't': [0, 1, 0],
        'f': [1, 0, 0]
    })
    tker = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .30, 'b': 2}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': 1}} 
    }
    lnd.updateTraps(traps, tker)

Doing this has the advantage (over creating a new landscape), of re-calculating only the parts of the landscape that are needed.
This can save significant amounts of computation in optimization routines.

To plot our updated landscape, we simply call our plotting routines again:

.. code-block:: python


    (fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
    lnd.plotSites(fig, ax[0])
    lnd.plotMigrationNetwork(fig, ax[0])
    lnd.plotTraps(fig, ax[0])
    lnd.plotTrapsNetwork(fig, ax[0])
    srv.plotMatrix(fig, ax[1], lnd.trapsMigration, lnd.trapsNumber)
    [srv.plotClean(fig, i, frame=False) for i in ax]


Where the traps with the black outline are immovable.


.. image:: ../../img/demo_updatedLandscape.jpg

The code used for this tutorial can be found `in this link <https://github.com/Chipdelmal/MGSurvE/blob/main/MGSurvE/demos/Demo_XY.py>`_
