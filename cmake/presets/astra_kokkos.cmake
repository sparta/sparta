#################### BEGIN SPARTA OPTIONS ####################
#get_property(SPARTA_CXX_COMPILE_FLAGS_DOCSTRING CACHE SPARTA_CXX_COMPILE_FLAGS PROPERTY DOCSTRING)
set(SPARTA_CXX_COMPILE_FLAGS 
    -fPIC
    CACHE 
    STRING 
    "Compiler flags use when building .o files for spa_*"
    )
set(SPARTA_CMAKE_CONFIG_STRING
    astra_kokkos
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    )
set(BUILD_MPI_TPL
    ON
    CACHE
    BOOL
    "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
    )
set(PKG_MPI_STUBS
    OFF
    CACHE
    BOOL
    "Enable or disable sparta mpi stubs package. Default: OFF.")
#################### END SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER "mpicxx" CACHE STRING "")
set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-g -O3" CACHE STRING "")
set(CMAKE_AR "ar" CACHE STRING "")
set(CMAKE_SHARED_LINKER_FLAGS "-fPIC -shared" CACHE STRING "")
#################### END CMAKE OPTIONS ####################

#################### BEGIN MPI OPTIONS ####################
set(CMAKE_CXX_FLAGS "")
#################### END   MPI OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ENABLE_OpenMP ON CACHE STRING "")
#################### END   KOKKOS OPTIONS ####################
