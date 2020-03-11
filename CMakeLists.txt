# TODO: Is 3.13 widley available??
cmake_minimum_required(VERSION 3.13.0)
project(SPARTA)

# TODO: process SpArTA variables?

## ENABLE variables
## ENABLE variables will trigger a find package. If the find package fails, it
## should suggest setting env vars to the user's local installation of using a
## module to set env vars.
## EXAMPLES:
# SPARTA_ENABLE_MPI:BOOL
# SPARTA_ENABLE_MPICH:BOOL
# SPARTA_ENABLE_OMPI:BOOL
# SPARTA_ENABLE_KOKKOS:BOOL
# SPARTA_ENABLE_MPICH:BOOL
# SPARTA_ENABLE_FFT:BOOL
# SPARTA_ENABLE_OPENMP:BOOL

# TODO: process enable variables and set __CMAKE_SPARTA_CONFIG string?
# TODO: Use TriBits? Would this help with TPL enables?

## Machine files
## Machine files can be specified via cmake -C machine.file, they will cache
## variables (as defaults) that can be optionally overriden. If a set comes
## after set-cache, the cache is ignored.
## EXAMPLES:
# chama.cmake
# cmake -C chama.cmake
# cmake -C chama.cmake -DSPARTA_ENABLE_MPICH
# cmake -C chama.cmake -DSPARTA_BIGBIG

add_subdirectory(src)
#add_subdirectory(lib)?
