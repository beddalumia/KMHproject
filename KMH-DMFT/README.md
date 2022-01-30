Dynamical-Mean-Field-Theory treatment of the Kane-Mele-Hubbard model.
Based on [QcmPlab](https://github.com/QcmPlab) libraries.

--------

- [KMH-DMFT_f90](./KMH-DMFT_f90) contains the main `Fortran` programs, together with an appropriate makefile.

- [KMH-DMFT_mat](./KMH-DMFT_mat) collects some `MATLAB` scripts, providing various workflows for interaction-driven transition lines. Based on [DMFT-LAB](https://github.com/bellomia/DMFT-LAB).

- [KMH-DMFT_f2py](./KMH-DMFT_f2py) will be the home for an interface to the `Python` package [PhaseMap](https://github.com/greschd/PhaseMap), allowing for smart & fast evaluations of multiparametric phase diagrams. Thorough testing of the python API provided in the [QcmPlab DMFT-ED library](https://github.com/QcmPlab/LIB_DMFT_ED) is a prerequisite. [**currently there is no stable API**]

- [KMH-DMFT_hpc](./KMH-DMFT_hpc) gives HPC workflows, interfacing all the scripts with the `SLURM` resource manager (intended for [Ulysses](https://www.itcs.sissa.it/services/computing/hpc), primarly).
