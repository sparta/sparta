if(SPARTA_DSMC_TESTING_PATH)
  get_filename_component(SPARTA_DSMC_TESTING_ABSOLUTE_PATH
                         ${SPARTA_DSMC_TESTING_PATH} ABSOLUTE)
  # message("dsmc_testing=${SPARTA_DSMC_TESTING_ABSOLUTE_PATH}")
  set(SPARTA_TEST_DRIVER python
                         ${SPARTA_DSMC_TESTING_ABSOLUTE_PATH}/regression.py)
endif()

# cmake-format: off
#
# sparta_add_test: Add a test called
# "${SPARTA_MACHINE}.${sparta_in_file}.mpi_${mpi_ranks}"
# @param sparta_in_file: The input file for the test being 
#                        added argument to sparta's "-in" option.
# @param mpi_ranks: A positive integer which is ignored if 
#                   BUILD_MPI is false.
#
# cmake-format: on
function(sparta_add_test sparta_in_file mpi_ranks)
  if(NOT BUILD_MPI)
    set(__test_name ${SPARTA_MACHINE}.${sparta_in_file})
    set(__sparta_command $<TARGET_FILE:${TARGET_SPARTA}>)
  else()
    set(__test_name ${SPARTA_MACHINE}.${sparta_in_file}.mpi_${mpi_ranks})
    set(__sparta_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${mpi_ranks}
                         $<TARGET_FILE:${TARGET_SPARTA}>)
  endif()

  if(SPARTA_DSMC_TESTING_PATH)
    # message("Adding test \"${__test_name}\" with test driver!")
    string(REPLACE ";" " " __sparta_driver_command "${__sparta_command}")
    set(__sparta_test_driver_postfix_args
        ${CMAKE_CURRENT_BINARY_DIR} -logread ${CMAKE_CURRENT_BINARY_DIR} olog
        -customtest ${sparta_in_file})
    add_test(
      NAME ${__test_name}
      COMMAND ${SPARTA_TEST_DRIVER} mpi_${mpi_ranks}
              "${__sparta_driver_command}" ${__sparta_test_driver_postfix_args})
    set_tests_properties(
      ${__test_name} PROPERTIES PASS_REGULAR_EXPRESSION "passed;no failures"
                                FAIL_REGULAR_EXPRESSION "FAILED")
  else()
    # message("Adding test \"${__test_name}\" without test driver!")
    add_test(NAME ${__test_name}
             COMMAND ${__sparta_command} -in ${sparta_in_file} -log
                     ${SPARTA_MACHINE}.${sparta_in_file}.log)
    set_tests_properties(
      ${__test_name}
      PROPERTIES PASS_REGULAR_EXPRESSION "" FAIL_REGULAR_EXPRESSION
                 "Error;ERROR;exited on signal")
  endif()
endfunction()

# cmake-format: off
#
# sparta_add_tests: Add the tests listed in in_file_list for all ranks listed in
# mpi_ranks
# @param in_file_list: a list of in.* files
# @param mpi_ranks:    a list of positive integers or "none"
#
# cmake-format: on
function(sparta_add_tests in_file_list mpi_ranks)
  foreach(mpi_rank IN LISTS mpi_ranks)
    foreach(in_file IN LISTS in_file_list)
      # message("sparta_add_test(${in_file} ${mpi_rank})")
      sparta_add_test(${in_file} ${mpi_rank})
    endforeach()
  endforeach()
endfunction()

# cmake-format: off
#
# sparta_add_all_tests: Add all the tests (*.in) in the current working
# directory with 1 and 4 mpi ranks
#
# cmake-format: on
function(sparta_add_all_tests mpi_ranks)
  file(
    GLOB __in_file_list
    LIST_DIRECTORIES false
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    CONFIGURE_DEPENDS in.*)

  if(BUILD_MPI)
    sparta_add_tests("${__in_file_list}" "${mpi_ranks}")
  else()
    sparta_add_tests("${__in_file_list}" "none")
  endif()
endfunction()
