# TODOs:
# set(CMAKE_OSX_SYSROOT "" CACHE STRING "")
#
# Seems like kitware doesn't want anyone using ar flags other than qc
# set(CMAKE_STATIC_LINKER_FLAGS "-rc" CACHE STRING "")
#
# md5 of executable doesn't match. Not clear if this is due to .a and exe file
# name differences between the old build and this cmake build. Need to verify.

#get_property(SPARTA_CXX_COMPILE_FLAGS_DOCSTRING CACHE SPARTA_CXX_COMPILE_FLAGS PROPERTY DOCSTRING)
set(SPARTA_CXX_COMPILE_FLAGS 
    -fPIC -DSPARTA_UNORDERED_MAP
    CACHE 
    STRING 
    "Compiler flags use when building .o files for spa_*"
    )

set(CMAKE_C_COMPILER "cc" CACHE STRING "")
set(CMAKE_CXX_COMPILER "c++" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-O" CACHE STRING "")

set(CMAKE_AR "ar" CACHE STRING "")

set(CMAKE_SHARED_LINKER_FLAGS "-fPIC -shared" CACHE STRING "")

set(Kokkos_ENABLE_OpenMP OFF CACHE STRING "")
