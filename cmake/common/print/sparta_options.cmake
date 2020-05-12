message(STATUS "")
message(STATUS "****************** Sparta Settings ******************")

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
