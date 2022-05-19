# Paper's Sims

This simulations set is generated as a companion for the initial publication of our framework.

## Requirements

* MGSurvE with dependencies [installed](https://chipdelmal.github.io/MGSurvE/build/html/installation.html).
* `./sims_out` folder in the path where the scripts are being run (for output).

## Running the scripts

Each landscape generator can be run independently for the desired demonstration in the following way:

```bash
python Landscape.py Ring
python Optimization.py Ring_LND_HOM
python Optimization.py Ring_LND_HET
```

just swapping the name of the landscape for one of the valid options ('Grid', 'Uniform', 'Ring', 'Poisson').

Alternatively, the whole set of experiments can be run with:

```bash
./PaperExperiments.sh
```

after giving it executable permissions:

```bash
chmod +x PaperExperiments.sh
```

## Running the scripts (Docker)

