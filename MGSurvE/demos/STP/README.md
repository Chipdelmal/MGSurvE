# São Tomé Demo


## [STP-Discrete](./STP-Discrete.py)

Runs the main optimization cycle. Inputs for the file are:

* `FXD_TRPS`: If `True`, two traps are set as immovable; otherwise, they are all optimizible.
* `AP`: Optimization summary stat identifier
    * `'man'` uses `np.mean`
    * `'sum'` uses `np.sum`
    * `'max'` uses `np.max`
* `TRPS_NUM`: Number of traps to distribute.
* `RID`: Run ID number (for repeated iterations of the cycle).

To re-run the whole dataset shown in our publication, run `bash STP_Discrete.sh` but beware that this is a computationally-intensive task that will likely take several days to finish.

## [STP-Compare](./STP-Compare.py)

Runs the analysis out of the optimization cycle. Inputs for the file follow:

* `ID`: Landscape identifier (`STPD` for São Tomé Discrete)
* `AP`: Optimization summary stat identifier
    * `'man'` uses `np.mean`
    * `'sum'` uses `np.sum`
    * `'max'` uses `np.max`


## [STP-Video](./STP-Video.py)

Generates the frames required to generate the videos that accompany our publication. Running `bash STP-Videos.sh` will export a video for each repetition of each optimization cycle, but is also quite computationally-expensive and will take several days to finish.