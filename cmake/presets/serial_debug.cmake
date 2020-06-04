# serial_debug = Linux box, g++, serial (no MPI)

include(${CMAKE_CURRENT_LIST_DIR}/serial.cmake)

#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    serial_debug
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g" CACHE STRING "" FORCE)
#################### END CMAKE OPTIONS ####################
