# chama_big = mpic++, OpenMPI (with default modules), BIGBIG setting

include(${CMAKE_SOURCE_DIR}/cmake/presets/chama.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    chama_big
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
set(SPARTA_CXX_COMPILE_FLAGS 
    ${SPARTA_CXX_COMPILE_FLAGS} -DSPARTA_BIGBIG
    CACHE 
    STRING 
    "Compiler flags use when building .o files for spa_*"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################
