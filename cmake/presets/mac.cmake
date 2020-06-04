# mac = Apple PowerBook laptop, c++, serial (no MPI)

# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    mac
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration")
set(SPARTA_CXX_COMPILE_FLAGS
    -fPIC -DSPARTA_UNORDERED_MAP
    CACHE STRING "Compiler flags use when building .o files for spa_*")
# ################### END SPARTA OPTIONS ####################

# ################### BEGIN TPL OPTIONS ####################
set(BUILD_MPI
    OFF
    CACHE
      BOOL
      "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
)
# ################### END   TPL OPTIONS ####################

# ################### BEGIN PACKAGE OPTIONS ####################
set(PKG_MPI_STUBS
    ON
    CACHE BOOL "Enable or disable sparta mpi stubs package. Default: OFF.")
# ################### END   PACKAGE OPTIONS ####################

# ################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER
    "cc"
    CACHE STRING "")
set(CMAKE_CXX_COMPILER
    "c++"
    CACHE STRING "")
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -O"
    CACHE STRING "")
set(CMAKE_AR
    "ar"
    CACHE STRING "")
set(CMAKE_SHARED_LINKER_FLAGS
    "-fPIC -shared"
    CACHE STRING "")
# ################### END CMAKE OPTIONS ####################

# ################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ENABLE_OPENMP
    OFF
    CACHE STRING "")
# ################### END   KOKKOS OPTIONS ####################
