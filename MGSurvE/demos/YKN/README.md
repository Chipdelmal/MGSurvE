# Yorkeys Knob Demo


## [YKN-Continuous](./YKN-Continuous.py) and [YKN-Discrete](./YKN-Discrete.py)

Runs the main optimization cycle. Inputs for the file are:

* `ID`: Landscape identifier (`YKN` for this landscape).
* `AP`: Optimization summary stat identifier:
    * `'man'` uses `np.mean`
    * `'sum'` uses `np.sum`
    * `'max'` uses `np.max`
* `RID`: Run ID number (for repeated iterations of the cycle).

To re-run the whole dataset shown in our publication, run `bash YKN_Discrete.sh` and `bash YKN_Continuous.sh` but beware that this is a computationally-intensive task that will likely take several days to finish.

## [YKN-Compare](./YKN-Compare.py)

Runs the analysis out of the optimization cycle. Inputs for the file follow:

* `ID`: Landscape identifier (`YKNC` for Yorkeys Knob Continuous and `YKND` for Yorkeys Knob Discrete)
* `AP`: Optimization summary stat identifier
    * `'man'` uses `np.mean`
    * `'sum'` uses `np.sum`
    * `'max'` uses `np.max`


## [YKN-Video](./YKN-Video.py)

Generates the frames required to generate the videos that accompany our publication. Running `bash YKN-Videos.sh` will export a video for each repetition of each optimization cycle, but is also quite computationally-expensive and will take several days to finish.