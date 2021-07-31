Dynamical-Mean-Field-Theory treatment of the Kane-Mele-Hubbard model.
Based on [QcmPlab](https://github.com/QcmPlab) libraries.

--------

- [KMH-DMFT_f90](./KMH-DMFT_f90) contains the main fortran programs, together with an appropriate makefile.

- [KMH-DMFT_mat](./KMH-DMFT_mat) defines a very basic set of workflows for exploring interaction-driven transitions, to different degrees of "automatization" for the control of simulation parameters; within `MATLAB`.

- [KMH-DMFT_py](./KMH-DMFT_py) will be the home for an interface to the `Python` package [PhaseMap](https://github.com/greschd/PhaseMap), allowing for smart&fast evaluations of multiparametric phase diagrams. Thorough testing of the python API provided in the [QcmPlab DMFT library](https://github.com/QcmPlab/LIB_DMFT_ED) is a prerequisite.

- [KMH-DMFT_jl](./KMH-DMFT_jl) will host a `julia-language` powered boost of the KMH-DMFT_py workflow. Hopefully f90-like performance.

- [KMH-DMFT_sh](./KMH-DMFT_sh) is a database of oldie-but-goldie `BASH` scripts, controlling traditional phase diagram spans & the like.

- [KMH-DMFT_hpc](./KMH-DMFT_hpc) gives HPC workflows, interfacing all the scripts with the `SLURM` resource manager (intended for [Ulysses](https://www.itcs.sissa.it/services/computing/hpc), primarly).
