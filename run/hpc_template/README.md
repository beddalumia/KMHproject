# KMH-HPC
Running QcmPlab DMFT codes on SLURM-managed hpc clusters

-------------

Everything is heavily based on the mighty [How to Survive the Mermaids](https://ulysses.readthedocs.io/index.html) unofficial guide to Ulysses (the SISSA hpc facility), with some slight twist based on the alerts reported in the [official-but-hyperconcise docs](https://www.itcs.sissa.it/services/computing/hpc), provided by ITCS.

- `secli.sh` is the original template given in the guide.
- `full-mpi_single-line_matjob.sh` is a preconfigured script for MPI jobs with variable number of nodes (default to 1). It calls the `+runDMFT` functions.
- `mpi-serial_single-line_matjob.sh` is a preconfigured script for MPI jobs with a _fixed_ mpi-rank (set to 1), to give _serial_ jobs. It calls the `+runDMFT` functions.
- `post_full-diagram_matjob.sh` is a preconfigured _serial_ job that calls the `PostDMFT.m` script, recollects all the `.mat` files, preserving the folder structure, and packs a tarball to be shipped by `scp`.




