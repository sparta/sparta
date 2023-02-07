# Compiling and installing sparta
```bash
cd /path/to/sparta

# Create a clean in-source build directory
rm -rf build
mkdir build install
cd build

# Configure the build system
cmake -DCMAKE_INSTALL_PREFIX=../install /path/to/sparta/cmake

# List project-specific options
cmake -L /path/to/sparta/cmake

# List project-specific options and help strings
cmake -LH /path/to/sparta/cmake

# List all options and help strings
cmake -LAH

# Show all generated targets
make help

# Build the project
make

# Install the binaries and libraries
make install
```
# SPARTA preset files
The sparta build system provides commonly used configuration settings as cmake files under
`/path/to/sparta/cmake/presets`. This allows developers to quickly set
configuration knobs before building sparta.
## Examples
### Using a preset file
```bash
cmake -C /path/to/<NAME>.cmake /path/to/sparta/cmake
```

# SPARTA Keyword Listing
## Package options
* PKG_FFT
  * Whether to enable the SPARTA FFT package.
* PKG_KOKKOS
  * Whether to enable the SPARTA KOKKOS package.
* PKG_MPI_STUBS
  * Whether to enable the SPARTA MPI_STUBS package.

## Third Party Library (TPL) options
* BUILD_KOKKOS
  * Whether to enable the SPARTA KOKKOS TPL.
* BUILD_MPI
  * Whether to enable the SPARTA MPI TPL.
* BUILD_JPEG
  * Whether to enable the SPARTA JPEG TPL.
* BUILD_PNG
  * Whether to enable the SPARTA PNG TPL.
* BUILD_MPI
  * Whether to enable the SPARTA MPI TPL.
* FFT
  * Which SPARTA FFT TPL to enable: FFTW2, FFTW3, or MKL.

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
  * Print the SPARTA packages and exit.
* SPARTA_LIST_TPLS
  * Print the SPARTA TPLs and exit.
* SPARTA_ENABLE_TESTING
  * Add tests in examples to be run via ctest.
* SPARTA_ENABLE_PARAVIEW_TESTING
  * Enable ParaView tests. Default if OFF.
  * When ON, must specify SPARTA_PARAVIEW_BIN_DIR and SPARTA_PARAVIEW_MPIEXEC.
* SPARTA_PARAVIEW_BIN_DIR
  * Path to ParaView bin directory containing pvbatch and pvpython.
* SPARTA_PARAVIEW_MPIEXEC
  * Path to program used to start ParaView mpi jobs, typically mpiexec in SPARTA_PARAVIEW_BIN_DIR.
* SPARTA_DSMC_TESTING_PATH
  * Add tests in SPARTA_DSMC_TESTING_PATH/examples to be run via ctest.
  * Run all tests via SPARTA_DSMC_TESTING_PATH/regression.py.
* SPARTA_SPA_ARGS
  * Additional arguments for the sparta binary. Only applied if SPARTA_ENABLE_TESTING or
  SPARTA_DSMC_TESTING_PATH are enabled.
* SPARTA_DSMC_TESTING_DRIVER_ARGS
  * Additional arguments for SPARTA_DSMC_TESTING_PATH/regression.py.
* SPARTA_CTEST_CONFIGS
  * Additional ctest configurations, separated by `;`, that allow `SPARTA_SPA_ARGS_<CONFIG_NAME>` or `SPARTA_DSMC_TESTING_DRIVER_ARGS_<CONFIG_NAME>` to be specified.
* SPARTA_MULTIBUILD_CONFIGS
  * Additional build configurations, separated by `;`, build sparta with the cache file from `SPARTA_MULTIBUILD_PRESET_DIR/<CONFIG_NAME>.cmake`.
* SPARTA_MULTIBUILD_PRESET_DIR
  * The path to custom preset files when using `SPARTA_MULTIBUILD_CONFIGS`. Only applied if `SPARTA_MULTIBUILD_CONFIGS` is enabled.

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

### Specifying multiple ctest configurations via the command line
```bash
cmake -DSPARTA_CTEST_CONFIGS="PARALLEL;SERIAL" \
      -DSPARTA_SPA_ARGS_SERIAL=spa_args \
      -DSPARTA_DSMC_TESTING_DRIVER_ARGS_PARALLEL=driver_args \
      /path/to/sparta/cmake

make -j

ctest -C SERIAL
ctest -C PARALLEL
```

### Specifying multiple build configurations via the command line
```bash
# Assumes that /path/to/sparta/cmake/presets/{test_mac_mpi,test_mac}.cmake exist
cmake -DSPARTA_MULTIBUILD_CONFIGS="test_mac;test_mac_mpi" \
      -DSPARTA_MULTIBUILD_PRESET_DIR=/path/to/sparta/cmake/presets/ \
      /path/to/sparta/cmake

make -j

ctest -VV
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
cmake --log-level=VERBOSE [-C /path/to/sparta/cmake/presets/<NAME>.cmake] /path/to/sparta/cmake
make VERBOSE=1
```
