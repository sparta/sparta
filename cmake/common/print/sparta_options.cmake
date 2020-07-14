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

if(FFT)
  message(STATUS "  ${FFT}")
endif()

message(STATUS "")

message(STATUS "Enabled extra options")
foreach(opt IN LISTS SPARTA_EXTRA_OPTIONS_LIST)
  if(${${opt}})
    message(STATUS "  ${opt}")
  endif()
endforeach()

message(STATUS "")

if(SPARTA_ENABLE_TESTING AND NOT SPARTA_DSMC_TESTING_PATH)
  message(STATUS "Test correctness verfication is limited without SPARTA_DSMC_TESTING_PATH set!")
endif()

message(STATUS "SPARTA_MACHINE: ${SPARTA_MACHINE}")
message(
  STATUS "SPARTA_DEFAULT_CXX_COMPILE_FLAGS: ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS}"
)
message(STATUS "SPARTA_CXX_COMPILE_FLAGS: ${SPARTA_CXX_COMPILE_FLAGS}")

# ######### BEGIN CHECK FOR OLD BUILD ARTIFACTS ##########
message(STATUS "Checking for old sparta build artifacts")
execute_process(COMMAND ${SPARTA_CMAKE_DIR}/check_for_sparta_build_artifacts.sh
                        ${SPARTA_SRC_DIR} OUTPUT_VARIABLE GNUMAKE_STDOUT)

if(DEFINED GNUMAKE_STDOUT)
  if(NOT ${GNUMAKE_STDOUT} STREQUAL "")
    message(FATAL_ERROR "${GNUMAKE_STDOUT}")
  endif()
endif()
message(STATUS "Checking for old sparta build artifacts - done")
# ######### END CHECK FOR OLD BUILD ARTIFACTS ##########

message(STATUS "****************** End   SPARTA Settings ******************")
message(STATUS "")
