GA Optimization
------------

In this demo, we will be optimizing the traps' positions to minimize the time it takes for a mosquito to get caught.
This is done with the `DEAP package <https://deap.readthedocs.io/en/master/>`_, as it allows much flexibility and implementation speedups.

We start by defining the parameters for our genetic algorithm:

.. code-block:: python

    (GENS, VERBOSE) = (2000, True)
    POP_SIZE = int(10*(lnd.trapsNumber*1.25))
    MAT = {'mate': .5, 'cxpb': 0.5}, 
    MUT = {'mean': 0, 'sd': max([i[1]-i[0] for i in bbox])/4, 'mutpb': .5, 'ipb': .5},
    SEL = {'tSize': 3}


Next, as defined by the `DEAP docs <https://deap.readthedocs.io/en/master/examples/index.html>`_, we register all the functions and operations
that we are going to use in our optimization cycle. For this version, we'll be using a pretty "vanilla" GA with
cxBlend, gaussian mutation, and tournament selection.

.. code-block:: python

    toolbox = base.Toolbox()
    creator.create("FitnessMin", base.Fitness, 
        weights=(-1.0, )
    )
    # Population creation -----------------------------------------------------
    creator.create(
        "Individual", list, 
        fitness=creator.FitnessMin
    )
    toolbox.register(
        "initChromosome", srv.initChromosome, 
        trapsCoords=lndGA.trapsCoords, 
        fixedTrapsMask=trpMsk, coordsRange=bbox
    )
    toolbox.register(
        "individualCreator", tools.initIterate, 
        creator.Individual, toolbox.initChromosome
    )
    toolbox.register(
        "populationCreator", tools.initRepeat, 
        list, toolbox.individualCreator
    )
    # Mutation and Crossover --------------------------------------------------
    toolbox.register(
        "mate", srv.cxBlend, 
        fixedTrapsMask=trpMsk, alpha=MAT['mate']
    )
    toolbox.register(
        "mutate", srv.mutateChromosome, 
        fixedTrapsMask=trpMsk, 
        randArgs={'loc': MUT['mean'], 'scale': MUT['sd']}
    )
    # Select and evaluate -----------------------------------------------------
    toolbox.register(
        "select", tools.selTournament, 
        tournsize=SEL['tSize']
    )
    toolbox.register(
        "evaluate", srv.calcFitness, 
        landscape=lndGA,
        optimFunction=srv.getDaysTillTrapped,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max}
    )

It is important to note that we provide custom implementations for the :code:`initChromosome`, :code:`cxBlend`, and :code:`mutateChromosome`; 
to allow immovable traps to be laid in the landscape, but we will stick to `DEAP's' <https://deap.readthedocs.io/en/master/>`_ implementations for this first exercise.

We now register summary statistics for our algorithm:

.. code-block:: python

    pop = toolbox.populationCreator(n=POP_SIZE)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)   
    stats.register("min", np.min)
    stats.register("avg", np.mean)
    stats.register("max", np.max)
    stats.register("traps", lambda fitnessValues: pop[fitnessValues.index(min(fitnessValues))])
    stats.register("best", lambda fitnessValues: fitnessValues.index(min(fitnessValues)))


Where the statistics go as follow (more stats can be added as needed):

* min: Traps' population minimum fitness (best in generation).
* avg: Traps' population average fitness.
* max: Traps' population maximum fitness (worst in generation).
* traps: Best traps positions in the current generation.
* best: Best fitness across populations.

Now, we run our optimization cycle:

.. code-block:: python

    (pop, logbook) = algorithms.eaSimple(
        pop, toolbox, cxpb=MAT['cxpb'], mutpb=MUT['mutpb'], ngen=GENS, 
        stats=stats, halloffame=hof, verbose=VERBOSE
    )

This will take some time depending on the number of generations and the size of the landscape/traps but once it's done running, we can get our resulting optimized positions by running:

.. code-block:: python

    bestChromosome = hof[0]
    bestPositions = np.reshape(bestChromosome, (-1, 2))