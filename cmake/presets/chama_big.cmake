include(${CMAKE_SOURCE_DIR}/cmake/presets/chama.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
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
