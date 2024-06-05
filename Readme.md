# RSOPT

 This is a package for optimization of S-wave velocity model based on surface wave dispersion and receiver function data. The package employs the non-linear iterative optimization algorithm to minimize the misfit between observed and predicted dispersion curves and receiver functions.

## Installation

### Install dependencies via conda

```bash
conda create -n rsopt -c conda-forge cxx-compiler fortran-compiler hdf5 cmake
conda activate rsopt
```

### Download the source code

```bash
git clone https://github.com/xumi1993/RSOpt.git
cd RSOpt
```

### Compile source code

```bash

mkdir build && cd build
cmake .. && make -j
```
