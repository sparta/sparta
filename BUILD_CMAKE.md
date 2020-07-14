# Compiling and installing sparta
```bash
cd /path/to/sparta

# Create a clean in-source build directory
rm -rf build
mkdir build install
cd build

# Configure the build system
cmake -DCMAKE_INSTALL_PREFIX=../install ../cmake

# List project-specific options
cmake -L ../cmake

# List project-specific options and help strings
cmake -LH ../cmake

# List all options and help strings
cmake -LAH

# Show all generated targets
make help

# Build the project
make

# Install the binaries and libraries
make install
```
# Sparta preset files
The sparta build system provides commonly used configuration settings as cmake files under
`/path/to/sparta/cmake/presets`. This allows developers to quickly set
configuration knobs before building sparta.
## Examples
### Using a preset file
```bash
cmake -C /path/to/<NAME>.cmake /path/to/sparta/cmake
```

# Sparta Keyword Listing
## Package options
* PKG_FFT
  * Whether to enable the Sparta FFT package.
* PKG_KOKKOS
  * Whether to enable the Sparta KOKKOS package.
* PKG_MPI_STUBS
  * Whether to enable the Sparta MPI_STUBS package.

## Third Party Library (TPL) options
* BUILD_KOKKOS
  * Whether to enable the Sparta KOKKOS TPL.
* BUILD_MPI
  * Whether to enable the Sparta MPI TPL.
* BUILD_JPEG
  * Whether to enable the Sparta JPEG TPL.
* BUILD_PNG
  * Whether to enable the Sparta PNG TPL.
* BUILD_MPI
  * Whether to enable the Sparta MPI TPL.
* FFT
  * Which Sparta FFT TPL to enable: FFTW2, FFTW3, or MKL.

Note: To point to a TPL installation, export <TPL>_ROOT=/path/to/tpl/install
before running cmake.

## Other options
* SPARTA_MACHINE
  * String to form the `spa_$SPARTA_MACHINE` binary file name.
* SPARTA_CXX_COMPILE_FLAGS
  * Selected compiler flags used when building object files for `spa_$SPARTA_MACHINE`.
* SPARTA_DEFAULT_CXX_COMPILE_FLAGS
  * Default compiler flags used when building object files for `spa_$SPARTA_MACHINE`.
* SPARTA_LIST_PKGS
  * Print the Sparta packages and exit.
* SPARTA_LIST_TPLS
  * Print the Sparta TPLs and exit.
* SPARTA_ENABLE_TESTING
  * Add tests in examples to be run via ctest.
* SPARTA_DSMC_TESTING_PATH
  * Add tests in SPARTA_DSMC_TESTING_PATH/examples to be run via ctest.
  * Run all tests via SPARTA_DSMC_TESTING_PATH/regression.py.
* SPARTA_SPA_ARGS
  * Additional arguments for the sparta binary. Only applied if SPARTA_ENABLE_TESTING or
  SPARTA_DSMC_TESTING_PATH are enabled.
* SPARTA_DSMC_TESTING_DRIVER_ARGS
  * Additional arguments for SPARTA_DSMC_TESTING_PATH/regression.py

## Examples
### Selecting packages via the command line
```bash
cmake -DPKG_<NAME>=[ON|OFF] /path/to/sparta/cmake
```

### Selecting TPLs via the command line
```bash
cmake -DBUILD_<NAME>_TPL=[ON|OFF] /path/to/sparta/cmake
```

### Specifying build flags via the command line
```bash
cmake -DSPARTA_DEFAULT_CXX_COMPILE_FLAGS=<FLAGS> /path/to/sparta/cmake
```

# Build system design
## Targets and dependency resolution
This build system consists of four targets:

1. `spa_$CONFIG_STRING`: The final sparta executable
2. `pkg_fft`: The optional FFT package
3. `pkg_mpi_stubs`: The optional MPI STUBS package
4. `pkg_kokkos`: The optional kokkos wrapper package

Every target is responsible for resolving its own dependencies. Every target A that
relies on another target B will pull in the dependencies that target B resolved.

Targets 2-4 are all optional packages that are built as static libraries.

Target 1 links against targets 2-4 (if enabled).

## The structure of the `sparta/cmake` directory
This directory contains two directories: `common` and `presets`.
### presets
Contains preset options that can be selected via: 
`cmake -C /path/to/presets/<NAME>.cmake`
### common
Contains three directories: `set`, `process`, and `print`.  Each of these
directories contains cmake files that are included by the top-level
`CMakeLists.txt` file in sparta. These `common` cmake files set build options,
process those build options, and finally print the settings that were selected.

# Build system triaging
## Quick start
```bash
cmake --log-level=VERBOSE [-C ../cmake/presets/<NAME>.cmake] ..
make VERBOSE=1
```
