# RunningDMFT
Running QcmPlab DMFT codes from MATLAB

----------

#### Choose your favorite workflow:

## Production

- `_dry.m` The most basic scenario: linear increase of U, no feedback, just some noninvasive check of convergence to flag out, _a posteriori_, which points should be discarded. 

- `_autostop.m` Basic linspace in U, without any feedback mechanism. Stops when dmft does not converge (to avoid wasting cpu-time).

- `_autostep.m` Gradually increases U by feedback-controlled steps, i.e. if dmft does not converge the point is discarded and the step reduced. Self-mixing is fixed (but you might update it manually, on the flight).


## Refinement

- `_refresh.m` Systematically enters the folder structure of a pre-existent calculation and uses the given restart and input files to perform some additional dmft loops. Useful if something has been added to the driver.

## Testing & Setup

- `_livemixing.m` Interactive-ish workflow: on-the-flight manual updates of the mixing ratio, while dmft waits for you (in a dumb way); Hubbard steps are inevitably fixed.


## Work-in-Progress & To-Do

- [ ] `_automixing.m` Automatic controll of self-mixing. Requires inspecting evolution of the dmft error, could be cumbersome.

- [ ] `_ai.m` Automatic control of both self-mixing and Hubbard steps. Difficult for sure, maybe just not worth all the effort.

- [ ] `_jobarray.m` Exploiting the array-env-variables provided by `SLURM` on HPC facilities.
