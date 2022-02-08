# KMHproject
A collection of programs and scripts to solve and analyze the Kane-Mele-Hubbard model in a variety of (dynamical-)mean-field settings. Based on the [QcmPlab](https://github.com/QcmPlab) libraries.

### Dependencies
* [SciFortran](https://github.com/QcmPlab/SciFortran.git), with all the dependencies listed therein. **[required]**
* [DMFT-tools](https://github.com/QcmPlab/DMFTtools.git), which depends on SciFortran. **[required]**
* [DMFT-ED](https://github.com/QcmPlab/LIB_DMFT_ED.git), which depends on SciFortran. **[required]**
* [DMFT-LAB](https://github.com/bellomia/DMFT-LAB), with a working installation of MATLAB / GNU Octave. **[required]**
* _[PhaseMap](https://github.com/greschd/PhaseMap), with all its [requirements](https://github.com/greschd/PhaseMap/blob/develop/setup.cfg)._[^1] 

All the **required** dependencies are [embedded](./Libraries/) in the project as [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules). To include them in the cloned repository run:

```
git clone --recursive https://github.com/bellomia/KMHproject.git KMHproj
```

To pull upstream changes of one specific submodule you can invoke:

```
git submodule update --remote --merge <submodule-name>
```

Or you can updated everything all together through the `foreach` syntax:

```
git submodule foreach `git pull origin`
git submodule update
```

**Warning:** The build and installation of the libraries is not automatized (yet?), thus you will need to enter each submodule directory and follow the provided instructions. All the upstream requirements have to be met, even the "optional" ones (e.g. MPI related stuff).

[^1]: Currently used only for the [KMH-MF_f2py](KMH-MF/KMH-MF_f2py) module. An alternative workflow, based on DMFT-LAB, is provided by the [KMH-MF_mat](KMH-MF/KMH-MF_mat) module.

