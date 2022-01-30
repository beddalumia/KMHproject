# KMHproject
A collection of programs and scripts to solve and analyze the Kane-Mele-Hubbard model in a variety of (dynamical-)mean-field settings. Based on the [QcmPlab](https://github.com/QcmPlab) libraries.

### Requirements
* [SciFortran](https://github.com/QcmPlab/SciFortran.git), with all the dependencies listed therein.
* [DMFT-tools](https://github.com/QcmPlab/DMFTtools.git), which depends on SciFortran.
* [DMFT-ED](https://github.com/QcmPlab/LIB_DMFT_ED.git), which depends on both SciFortran and DMFT-tools.
* [DMFT-LAB](https://github.com/bellomia/DMFT-LAB), with a working installation of MATLAB / GNU Octave.
* \< optional[^1] \> [PhaseMap](https://github.com/greschd/PhaseMap), with all its [requirements](https://github.com/greschd/PhaseMap/blob/develop/setup.cfg).


[^1]: Currently used only for the [KMH-MF_f2py](KMH-MF/KMH-MF_f2py) module. An alternative workflow, based on DMFT-LAB, is provided by the [KMH-MF_mat](KMH-MF/KMH-MF_mat) module.

