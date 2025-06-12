## Overview

This repository provides **direct evaluation of the free-space Stokeslet**, extracted from [SE_unified](https://github.com/joarbagge/SE_unified_v2) by Joar Bagge.  
It is a stripped-down version tailored specifically for the **free-space Stokes kernel**, without periodicity or other features.

## Build Instructions

A precompiled binary is included, tested with MATLAB R2025a on linux. To build from source **CMake** is used for cross-platform compilation.

To build the project:

```bash
cd build
cmake ..
make 
```

If CMake cannot find your Matlab installation, try

```
cmake .. -DMatlab_ROOT_DIR="/path/to/MATLAB/R20XXv"
```

with the path replaced by the path to the Matlab installation.
