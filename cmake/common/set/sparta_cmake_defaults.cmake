################################################################################
# This file sets common default options that all sparta builds use. 
# These options can be overridden at configure time via `cmake -DVAR=VAL` or
# `cmake -C /path/to/preset/presets.cmake`
################################################################################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS
    -DSPARTA_GZIP
    CACHE
    STRING "Compiler flags used when building object files for the \"spa_\" executable"
    )

set(SPARTA_MACHINE
    ""
    CACHE
    STRING "Suffix to append to spa binary (WON'T enable any features automatically)"
    )
