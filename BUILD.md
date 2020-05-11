# Overview
This build system relies on cmake targets to resolve dependencies. For SPARTA,
we have three types of targets:
    1. The final sparta executable
    2. Sparta packages: FFT, STUBS, KOKKOS
    3. Sparta TPLS: KOKKOS

Every target is responsible for knowing its dependencies. For example, the
SPARTA package "KOKKOS" depends on the TPL "KOKKOS" so the CMakeLists.txt file
for the SPARTA PACKAGE "KOKKOS" resolves this dependency.

## Quick Start
```bash
cd /path/to/sparta
mkdir build
cd build
cmake ..
cmake -LH
cmake -LAH
make
```

## Quick start build triaging
``bash
cmake --log-level=VERBOSE [-C ../cmake/presets/FILE.cmake] ..
make VERBOSE=1
```
