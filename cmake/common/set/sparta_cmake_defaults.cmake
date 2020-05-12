################################################################################
# This file sets common default options that all sparta builds use. 
# These options can be overridden at configure time via `cmake -DVAR=VAL` or
# `cmake -C /path/to/preset/presets.cmake`
################################################################################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS
    -DSPARTA_GZIP
    CACHE
    STRING "Compiler flags use when building .o files for spa_*"
    )

set(SPARTA_CMAKE_CONFIG_STRING
    nopreset
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    )
