# Yorkeys Knob Demo

Code required to replicate the "Discrete and continuous optimization on a suburban landscape (*Ae. aegypti* in Queensland, Australia)" demo of our publication. Please note that this example requires intermediate to advanced understanding of the MGSurvE framework. To get started with these concepts please have a look at our [documentation's tutorials](https://chipdelmal.github.io/MGSurvE/build/html/demos.html).

## Quick Instructions

To run these scripts download the whole folder (including the `GEO` subfolder, which contains the points dataset) to your local computer. After installing MGSurvE (including its [cartopy](https://scitools.org.uk/cartopy/docs/latest/) dependency), run the main file `YKN-Continuous.py`. This can be done either in an interactive session, or from the terminal as follows:

```bash
python YKN-Continuous.py 'YKN' 'max' 01
```

where the arguments correspond to the ones described in the following section of this document. The output files will be exported to a `sims_out` folder in the same directory by default.

## Files Description

### [YKN-Continuous](./YKN-Continuous.py) and [YKN-Discrete](./YKN-Discrete.py)

Runs the main optimization cycle. Inputs for the file are:

* `ID`: Landscape identifier (`YKN` for this landscape).
* `AP`: Optimization summary stat identifier:
    * `'man'` uses `np.mean`
    * `'sum'` uses `np.sum`
    * `'max'` uses `np.max`
* `RID`: Run ID number (for repeated iterations of the cycle).

To re-run the whole dataset shown in our publication, run `bash YKN_Discrete.sh` and `bash YKN_Continuous.sh` but beware that this is a computationally-intensive task that will likely take several days to finish.

### [YKN-Compare](./YKN-Compare.py)

Runs the analysis out of the optimization cycle. Inputs for the file follow:

* `ID`: Landscape identifier (`YKNC` for Yorkeys Knob Continuous and `YKND` for Yorkeys Knob Discrete)
* `AP`: Optimization summary stat identifier
    * `'man'` uses `np.mean`
    * `'sum'` uses `np.sum`
    * `'max'` uses `np.max`


### [YKN-Video](./YKN-Video.py)

Generates the frames required to generate the videos that accompany our publication. Running `bash YKN-Videos.sh` will export a video for each repetition of each optimization cycle, but is also quite computationally-expensive and will take several days to finish.


### [auxiliary](./auxiliary.py)

Contains functions that are not part of MGSurvE but which are required by the other scripts to run correctly.