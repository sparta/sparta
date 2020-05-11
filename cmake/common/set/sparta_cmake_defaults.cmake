################################################################################
# This file sets common default options that all sparta builds use. 
# These options can be overridden at configure time via `cmake -DVAR=VAL`.
################################################################################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS
    -DSPARTA_GZIP
    CACHE
    STRING "Compiler flags use when building .o files for spa_*"
    )
