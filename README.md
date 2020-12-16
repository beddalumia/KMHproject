Running [QcmPlab](https://github.com/QcmPlab) codes: a collection of scripts.

--------

- [RunningDMFT.m](./RunningDMFT.m) defines a very basic set of workflows for exploring interaction-driven Mott transitions, to different degrees of "automatization" for the control of simulation parameters; within `MATLAB`.

- [RunningDMFT.py](./RunningDMFT.py) will be the home for an interface to the `Python` package [PhaseMap](https://github.com/greschd/PhaseMap), allowing for smart&fast evaluations of multiparametric phase diagrams. Thorough testing of the python API provided in the [QcmPlab DMFT library](https://github.com/QcmPlab/LIB_DMFT_ED) is a prerequisite.

- [RunningDMFT.jl](./RunningDMFT.jl) will host a `julia-language` powered boost of the RunningDMFT.py workflow. Hopefully f90-like performance.

- [RunningDMFT.sh](./RunningDMFT.sh) is a database of oldie-but-goldie `BASH` scripts, controlling traditional phase diagram spans & the like.

- [RunningDMFT.slurm](./RunningDMFT.slurm) gives HPC workflows, interfacing all the scripts with the `SLURM` resource manager (intended for [Ulysses](https://www.itcs.sissa.it/services/computing/hpc), primarly).
