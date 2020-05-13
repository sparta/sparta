# serial = Linux box, g++, serial (no MPI)

#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    serial
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    )
set(SPARTA_CXX_COMPILE_FLAGS 
    -fPIC
    CACHE 
    STRING 
    "Compiler flags use when building .o files for spa_*"
    )
#################### END SPARTA OPTIONS ####################

#################### BEGIN TPL OPTIONS ####################
set(BUILD_MPI_TPL
    OFF
    CACHE
    BOOL
    "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
    FORCE
    )
#################### END   TPL OPTIONS ####################

#################### BEGIN PACKAGE OPTIONS ####################
set(PKG_MPI_STUBS
    ON
    CACHE
    BOOL
    "Enable or disable sparta mpi stubs package. Default: OFF."
    FORCE
    )
#################### END   PACKAGE OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER "gcc" CACHE STRING "")
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-fPIC -O3" CACHE STRING "")
set(CMAKE_AR "ar" CACHE STRING "")
set(CMAKE_SHARED_LINKER_FLAGS "-fPIC -shared" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS "-O" CACHE STRING "")
#################### END CMAKE OPTIONS ####################
