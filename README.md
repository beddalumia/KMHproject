# KMHproject
A collection of programs and scripts to solve and analyze the Kane-Mele-Hubbard model in a variety of (dynamical-)mean-field settings. Based on the [QcmPlab](https://github.com/QcmPlab) libraries.

### Dependencies

The project relies on several external libraries:

* [SciFortran](https://github.com/QcmPlab/SciFortran.git), with all the dependencies listed therein. **[required]**
* [DMFT-tools](https://github.com/QcmPlab/DMFTtools.git), which depends on SciFortran. **[required]**
* [DMFT-ED](https://github.com/QcmPlab/LIB_DMFT_ED.git), which depends on SciFortran. **[required]**
* [DMFT-LAB](https://github.com/bellomia/DMFT-LAB), requiring MATLAB or GNU Octave. **[required]**
* [PythTB](https://github.com/bellomia/PythTB), requiring Python > 2.7, numpy and matplotlib. **[optional]**

All the **required**[^1] dependencies are [embedded](./lib/) in the project as [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules). To download the project with all the intended versions of such libraries run:

```
git clone --recursive https://github.com/bellomia/KMHproject.git
```

To pull upstream changes < of one specific submodule > you can invoke:

```
git submodule update --remote --merge <submodule-name>
```

**Warning:** The build and installation of the libraries is not automatized (yet?), thus you will need to enter each submodule directory and follow the provided instructions. All the upstream requirements have to be met, even the "optional" ones (e.g. MPI related stuff).

[^1]: `PythTB` is used only as an authoritative benchmark for the noninteracting Kane Mele model, through the [`tb_km_2d.py`](src/tb_km_2d.py) script. The last stable version can be installed by typing `pip install pythtb --upgrade`.


