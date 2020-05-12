message(STATUS "")
message(STATUS "****************** Begin SPARTA Settings ******************")

message(STATUS "Enabled packages")
foreach(opt IN LISTS SPARTA_PKG_LIST)
    message(VERBOSE "opt=${${opt}}")
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
message(STATUS "****************** End   SPARTA Settings ******************")
