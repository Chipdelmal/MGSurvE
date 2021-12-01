GA Optimization
------------


.. code-block:: python

    (GENS, VERBOSE) = (2000, True)
    POP_SIZE = int(10*(lnd.trapsNumber*1.25))
    MAT = {'mate': .5, 'cxpb': 0.5}, 
    MUT = {'mean': 0, 'sd': max([i[1]-i[0] for i in bbox])/4, 'mutpb': .2, 'ipb': .2},
    SEL = {'tSize': 3}


.. code-block:: python

    lndGA = deepcopy(lnd)


.. code-block:: python
    
    toolbox = base.Toolbox()
    creator.create("FitnessMin", 
        base.Fitness, weights=(-1.0, )
    )
    # Population creation -----------------------------------------------------
    creator.create("Individual", 
        list, fitness=creator.FitnessMin
    )
    toolbox.register("initChromosome", srv.initChromosome, 
        trapsCoords=lndGA.trapsCoords, 
        fixedTrapsMask=trpMsk, coordsRange=bbox
    )
    toolbox.register("individualCreator", tools.initIterate, 
        creator.Individual, toolbox.initChromosome
    )
    toolbox.register("populationCreator", tools.initRepeat, 
        list, toolbox.individualCreator
    )
    # Mutation and Crossover --------------------------------------------------
    toolbox.register("mate", srv.cxBlend, 
        fixedTrapsMask=trpMsk, alpha=MAT['mate']
    )
    toolbox.register("mutate", srv.mutateChromosome, 
        fixedTrapsMask=trpMsk, 
        randArgs={'loc': MUT['mean'], 'scale': MUT['sd']}
    )
    # Select and evaluate -----------------------------------------------------
    toolbox.register("select", 
        tools.selTournament, tournsize=SEL['tSize']
    )
    toolbox.register("evaluate", 
        srv.calcFitness, 
        landscape=lndGA,
        optimFunction=srv.getDaysTillTrapped,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max}
    )