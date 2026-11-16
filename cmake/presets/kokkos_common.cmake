# kokkos_common = KOKKOS package with default settings

# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    kokkos_common
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration")
set(SPARTA_CXX_COMPILE_FLAGS
    -fPIC
    CACHE STRING "Compiler flags use when building .o files for spa_*")
# ################### END   SPARTA OPTIONS ####################

# ################### BEGIN TPL OPTIONS ####################
set(BUILD_MPI
    ON
    CACHE
      BOOL
      "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
)
# ################### END   TPL OPTIONS ####################

# ################### BEGIN PACKAGE OPTIONS ####################
set(PKG_MPI_STUBS
    OFF
    CACHE BOOL "Enable or disable sparta mpi stubs package. Default: OFF.")
set(PKG_KOKKOS
    ON
    CACHE BOOL "Enable or disable sparta kokkos package. Default: OFF.")
# ################### END   PACKAGE OPTIONS ####################

# ################### BEGIN CMAKE OPTIONS ####################
set(CRAYPE_VERSION $ENV{CRAYPE_VERSION})
if(CRAYPE_VERSION)
  set(CMAKE_C_COMPILER
      "cc"
      CACHE STRING "")
  set(CMAKE_CXX_COMPILER
      "CC"
      CACHE STRING "")
else()
  set(CMAKE_C_COMPILER
      "mpicc"
      CACHE STRING "")
  set(CMAKE_CXX_COMPILER
      "mpicxx"
      CACHE STRING "")
endif()
# optimization and debug flags come with the build type, as in LAMMPS: -g only
# for RelWithDebInfo (the default here) and Debug builds, not for Release
set(CMAKE_BUILD_TYPE
    RelWithDebInfo
    CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "-g -O3 -DNDEBUG"
    CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE
    "-O3 -DNDEBUG"
    CACHE STRING "")
set(CMAKE_AR
    "ar"
    CACHE STRING "")
set(CMAKE_SHARED_LINKER_FLAGS
    "-fPIC -shared"
    CACHE STRING "")
# ################### END CMAKE OPTIONS ####################

# ################### BEGIN MPI OPTIONS ####################
set(SPARTA_CXX_COMPILE_FLAGS
    ${SPARTA_CXX_COMPILE_FLAGS} -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX=1
    CACHE STRING "" FORCE)
# ################### END   MPI OPTIONS ####################

# ################### BEGIN KOKKOS OPTIONS ####################
# ################### END   KOKKOS OPTIONS ####################
