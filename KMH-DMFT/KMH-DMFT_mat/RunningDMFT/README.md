# KMH-DMFT_mat
Running QcmPlab DMFT codes from MATLAB

----------

#### Choose your favorite workflow

- `_dry.m` The most basic scenario: linear increase of U, no feedback, just some noninvasive check of convergence to flag out, _a posteriori_, which points should be discarded. 

- `_autostop.m` Basic linspace in U, without any feedback mechanism. Stops when dmft does not converge (to avoid wasting time).

- `_autoupdate.m` Gradually increases U by feedback-controlled steps, i.e. if dmft does not converge the point is discarded and the step reduced. Self-mixing is fixed (but you might update it manually, with some swag).

- `_livemixing.m` Interactive-ish workflow: on-the-flight manual updates of the mixing ratio, while dmft waits for you (in a dumb way); Hubbard steps are inevitably fixed.


#### To do

- [ ] `_automixing.m` Automatic controll of self-mixing. Requires inspecting evolution of the dmft error, could be cumbersome.

- [ ] `_ai.m` Automatic control of both self-mixing and Hubbard steps. Difficult for sure, maybe just not worth all the effort.

- [ ] `_SlurmArrays.m` Exploiting the array-env-variables provided by `SLURM` on HPC facilities. Totally broken at the moment, don't use it.
