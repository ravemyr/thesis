# The Spectral Ewald "Reset" package

## Building

The code is built with CMake. To build, open a terminal in the
root directory and do

```
cd build
cmake ..
make
```

If CMake cannot find your Matlab installation, try

```
cmake .. -DMatlab_ROOT_DIR="/path/to/MATLAB/R2019a"
```

with the path replaced by the path to the Matlab installation.
