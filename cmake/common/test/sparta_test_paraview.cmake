
set(SPARTA_PARAVIEW_BIN_DIR "" CACHE PATH "ParaView bin directory path")

set(SPARTA_PARAVIEW_MPIEXEC "" CACHE FILEPATH
    "Path to ParaView mpiexec command, if not in SPARTA_PARAVIEW_BIN_DIR")

set(PARAVIEW_TEST_NAME_REGEX ".paraview$")

function(is_paraview_test RETURN_VALUE TEST_NAME)
    if(${TEST_NAME} MATCHES ${PARAVIEW_TEST_NAME_REGEX})
        set(${RETURN_VALUE} TRUE PARENT_SCOPE)
    else()
        set(${RETURN_VALUE} FALSE PARENT_SCOPE)
    endif()
endfunction()

set(GRID2PARAVIEW_MODULE ${PROJECT_BINARY_DIR}/paraview/grid2paraview.py)
set(SURF2PARAVIEW_MODULE ${PROJECT_BINARY_DIR}/paraview/surf2paraview.py)

if(SPARTA_ENABLE_TESTING AND SPARTA_ENABLE_PARAVIEW_TESTING)

    if(SPARTA_PARAVIEW_BIN_DIR STREQUAL "")
        message(FATAL_ERROR
            "Must provide path to ParaView bin directory in SPARTA_PARAVIEW_BIN_DIR")
    endif()

    unset(PVPYTHON_EXECUTABLE CACHE)
    find_program(PVPYTHON_EXECUTABLE pvpython PATHS
        ${SPARTA_PARAVIEW_BIN_DIR} NO_DEFAULT_PATH)

    if(NOT PVPYTHON_EXECUTABLE)
        message(FATAL_ERROR
            "pvpython executable not found at path ${SPARTA_PARAVIEW_BIN_DIR}")
    endif()

    unset(PVBATCH_EXECUTABLE CACHE)
    find_program(PVBATCH_EXECUTABLE pvbatch PATHS
        ${SPARTA_PARAVIEW_BIN_DIR} NO_DEFAULT_PATH)

    if(NOT PVBATCH_EXECUTABLE)
        message(FATAL_ERROR
            "pvbatch executable not found at path ${SPARTA_PARAVIEW_BIN_DIR}")
    endif()

    if(SPARTA_PARAVIEW_MPIEXEC STREQUAL "")
        find_program(PARAVIEW_MPIEXEC mpiexec PATHS
            ${SPARTA_PARAVIEW_BIN_DIR} NO_DEFAULT_PATH)
        if(NOT PARAVIEW_MPIEXEC)
            message(FATAL_ERROR
                "mpiexec executable not found at path ${SPARTA_PARAVIEW_BIN_DIR}")
        endif()
        set(SPARTA_PARAVIEW_MPIEXEC ${PARAVIEW_MPIEXEC}
            CACHE PATH "ParaView mpiexec executable path" FORCE)
    endif()

endif()
