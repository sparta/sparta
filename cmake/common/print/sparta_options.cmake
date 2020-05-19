message(STATUS "")
message(STATUS "****************** Begin SPARTA Settings ******************")

message(STATUS "Enabled packages")
foreach(opt IN LISTS SPARTA_PKG_LIST)
    if(${${opt}})
        message(STATUS "  ${opt}")
    endif()
endforeach()

message(STATUS "")

message(STATUS "Enabled TPLs")
foreach(opt IN LISTS SPARTA_BUILD_TPL_LIST)
    if(${${opt}})
        message(STATUS "  ${opt}")
    endif()
endforeach()

message(STATUS "")

message(STATUS "SPARTA_DEFAULT_CXX_COMPILE_FLAGS: ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS}")
message(STATUS "SPARTA_CXX_COMPILE_FLAGS: ${SPARTA_CXX_COMPILE_FLAGS}")

########## BEGIN CHECK FOR OLD BUILD ARTIFACTS ##########
message(STATUS "Checking for old sparta build artifacts")
execute_process(
    COMMAND ${SPARTA_CMAKE_DIR}/check_for_sparta_build_artifacts.sh ${SPARTA_SRC_DIR}
    OUTPUT_VARIABLE GNUMAKE_STDOUT
)

if (DEFINED GNUMAKE_STDOUT)
    if (NOT ${GNUMAKE_STDOUT} STREQUAL "")
        message(FATAL_ERROR "${GNUMAKE_STDOUT}")
    endif()
endif()
message(STATUS "Checking for old sparta build artifacts - done")
########## END CHECK FOR OLD BUILD ARTIFACTS ##########

message(STATUS "****************** End   SPARTA Settings ******************")
