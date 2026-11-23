# openmpi = Linux box, g++, OpenMPI via explicit path

# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    openmpi
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration")
set(SPARTA_CXX_COMPILE_FLAGS
    -fPIC -DSPARTA_JPEG -ljpeg
    CACHE STRING "Compiler flags use when building .o files for spa_*")
# ################### END SPARTA OPTIONS ####################

# ################### BEGIN TPL OPTIONS ####################
# ################### END   TPL OPTIONS ####################

# ################### BEGIN PACKAGE OPTIONS ####################
# ################### END   PACKAGE OPTIONS ####################

# ################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER
    "gcc"
    CACHE STRING "")
set(CMAKE_CXX_COMPILER
    "g++"
    CACHE STRING "")
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
set(CMAKE_EXE_LINKER_FLAGS
    "-g -O"
    CACHE STRING "")
# ################### END CMAKE OPTIONS ####################

# ################### BEGIN MPI OPTIONS ####################
set(MPI_DIR "/usr/local/openmpi")
# ################### END MPI OPTIONS ####################
