# 🚧 Work ⚠️ in 🪜 Progress 🚧

- `KMH-MF-f2py` provides a tentative interface to [PhaseMap](https://github.com/greschd/PhaseMap), aiming at accurate-but-fast phase diagrams, with a number of evaluations that scales with the phase boundaries, instead of the full parameter space.    
The hartree-fock case is pretty straightforward and we have a working first implementation, but the DMFT poses a bunch of additional problems, especially:
    - the lack of a stable[^1] python API in [DMFT-ED](https://github.com/QcmPlab/LIB_DMFT_ED);
    - the need for a careful _restarting strategy_, in order to achieve easy convergence, which probably rules out the possibility of huge jumps in parameter space, which are instead somehow at the heart of PhaseMap's algorithm.

[^1]: [EDIpack](https://github.com/QcmPlab/EDIpack) has indeed a working API, but it solves only the normal phase in the _S<sub>z</sub>_ basis, so it is not suitable for the in-plane-magnetic phase of the KMH model.

- `PlottingTest.m` is a very rough script for phase-diagram plotting, whose core ideas should be merged soon into `+plotDMFT` module in [DMFT-LAB](https://github.com/bellomia/DMFT-LAB).
