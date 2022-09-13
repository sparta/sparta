find_package(Python REQUIRED)
if(SPARTA_DSMC_TESTING_PATH)
  get_filename_component(SPARTA_DSMC_TESTING_ABSOLUTE_PATH
                         ${SPARTA_DSMC_TESTING_PATH} ABSOLUTE)
  # message("dsmc_testing=${SPARTA_DSMC_TESTING_ABSOLUTE_PATH}")
  set(SPARTA_TEST_DRIVER Python::Interpreter
                         ${SPARTA_DSMC_TESTING_ABSOLUTE_PATH}/regression.py)
  set(PARAVIEW_TEST_DRIVER ${PVPYTHON_EXECUTABLE}
    ${SPARTA_DSMC_TESTING_ABSOLUTE_PATH}/paraview_regression.py)
else()
  set(SPARTA_TEST_DRIVER Python::Interpreter ${SPARTA_TOOLS_DIR}/testing/regression.py)
  set(PARAVIEW_TEST_DRIVER ${PVPYTHON_EXECUTABLE}
    ${SPARTA_TOOLS_DIR}/testing/paraview_regression.py)
endif()

# cmake-format: off
#
# sparta_add_test: Add a test called
# "${SPARTA_MACHINE}.${sparta_in_file}.mpi_${mpi_ranks}"
# @param sparta_in_file: The input file for the test being
#                        added argument to sparta's "-in" option.
# @param mpi_ranks:      A positive integer which is ignored if
#                        BUILD_MPI is false.
# @param config_name:    A string describing the configuration to add
#                        this test to. Setting config_name to "" will
#                        add the test to the default configuration.
# @param work_dir:       A optional working directory argument where
#                        the test will be run.
#
# cmake-format: on
function(sparta_add_test sparta_in_file mpi_ranks config_name)
  set(__total_params 4)
  if(${ARGC} EQUAL ${__total_params})
    if (NOT SPARTA_DSMC_TESTING_PATH)
      math(EXPR ARGC "${ARGC}-1")
      set(__work_dir "${ARGV${ARGC}}")
      get_filename_component(__work_dir_basename "${__work_dir}" NAME)
      set(__run_dir ${CMAKE_CURRENT_BINARY_DIR}/${__work_dir_basename})
    else()
      set(__run_dir ${CMAKE_CURRENT_BINARY_DIR})
    endif()
  else()
    set(__work_dir ${CMAKE_CURRENT_SOURCE_DIR})
    set(__run_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  # message("sparta_add_test: __work_dir=${__work_dir}")

  if(config_name STREQUAL "")
    set(__config_name "${config_name}")
  else()
    set(__config_name "_${config_name}")
    # The user must set SPARTA_SPA_ARGS_${config_name}
  endif()

  if(SPARTA_SPA_ARGS${__config_name})
    string(REPLACE " " ";" __spa_args ${SPARTA_SPA_ARGS${__config_name}})
  endif()
  if(NOT BUILD_MPI)
    set(__test_name ${SPARTA_MACHINE}.${sparta_in_file}${__config_name})
    set(__sparta_command $<TARGET_FILE:${TARGET_SPARTA}> ${__spa_args})
  else()
    set(__test_name
        ${SPARTA_MACHINE}.${sparta_in_file}.mpi_${mpi_ranks}${__config_name})
    set(__sparta_command
        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${mpi_ranks} --bind-to none
        $<TARGET_FILE:${TARGET_SPARTA}> ${__spa_args})
  endif()

  # message("Adding test \"${__test_name}\" with test driver!")
  string(REPLACE ";" " " __sparta_driver_command "${__sparta_command}")
  if(SPARTA_DSMC_TESTING_DRIVER_ARGS${__config_name})
    string(REPLACE " " ";" __driver_args
                   ${SPARTA_DSMC_TESTING_DRIVER_ARGS${__config_name}})
  endif()
  set(__sparta_test_driver_postfix_args
      ${__run_dir}
      -logread
      ${__run_dir}
      olog
      -customtest
      ${sparta_in_file}
      ${__driver_args})
  is_paraview_test(paraview_test sparta_in_file)
  if(paraview_test)
    add_test(
      NAME ${__test_name}
      CONFIGURATIONS ${config_name}
      COMMAND ${PARAVIEW_TEST_DRIVER} ${mpi_ranks} "${__sparta_driver_command}"
          ${sparta_in_file} ${SPARTA_PARAVIEW_MPIEXEC} ${PVPYTHON_EXECUTABLE}
          ${PVBATCH_EXECUTABLE} ${GRID2PARAVIEW_MODULE} ${SURF2PARAVIEW_MODULE}
          ${GRID2PARAVIEWCELLS_MODULE}
      WORKING_DIRECTORY ${__run_dir})
  else()
    add_test(
      NAME ${__test_name}
      CONFIGURATIONS ${config_name}
      COMMAND ${SPARTA_TEST_DRIVER} mpi_${mpi_ranks} "${__sparta_driver_command}"
              ${__sparta_test_driver_postfix_args}
      WORKING_DIRECTORY ${__run_dir})
    # unable to compile regex: ^\*{3} test .* passed
  endif()

  set_tests_properties(
    ${__test_name} PROPERTIES PASS_REGULAR_EXPRESSION "passed;no failures"
                              FAIL_REGULAR_EXPRESSION "FAILED")

  if(NOT SPARTA_DSMC_TESTING_THREADS_PER_RANK)
    # message(WARNING "SPARTA_DSMC_TESTING_THREADS_PER_RANK is uset! Defaulting
    # to 1.")
    set(SPARTA_DSMC_TESTING_THREADS_PER_RANK 1)
  endif()

  if(BUILD_MPI)
    math(EXPR processors
         "${mpi_ranks} * ${SPARTA_DSMC_TESTING_THREADS_PER_RANK}")
  else()
    set(processors ${SPARTA_DSMC_TESTING_THREADS_PER_RANK})
  endif()

  # message("processors=${processors}")
  set_property(TEST ${__test_name} PROPERTY PROCESSORS "${processors}")
endfunction()

# cmake-format: off
#
# sparta_add_tests_to_config: Add the tests listed in in_file_list for all ranks listed in
#                             mpi_ranks to the specified configuration.
# @param in_file_list: A list of in.* files
# @param mpi_ranks:    A list of positive integers or "none"
# @param config_name:  A string describing the configuration to add
#                      this test to. Setting config_name to "" will
#                      add the test to the default configuration.
# @param work_dir:     A optional working directory argument where
#                      the test will be run.
#
# cmake-format: on
function(sparta_add_tests_to_config in_file_list mpi_ranks config_name)
  # Too bad cmake uses ARGC/ARGV for macros too.
  set(__total_params 4)
  if(${ARGC} EQUAL ${__total_params})
    math(EXPR ARGC "${ARGC}-1")
    set(__work_dir "${ARGV${ARGC}}")
  else()
    set(__work_dir ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  # message("sparta_add_tests_to_config: __work_dir=${__work_dir}")

  foreach(mpi_rank IN LISTS mpi_ranks)
    foreach(in_file IN LISTS in_file_list)
      # message("sparta_add_test(${in_file} ${mpi_rank})")
      sparta_add_test(${in_file} ${mpi_rank} "${config_name}" ${__work_dir})
    endforeach()
  endforeach()
endfunction()

# Wrapper to sparta_add_tests_to_config
function(sparta_add_tests in_file_list mpi_ranks)
  sparta_add_test_to_config("${__in_file_list}" "${mpi_ranks}" "")
endfunction()

# cmake-format: off
#
# sparta_add_all_tests: Add all the tests (*.in) in the current working
#                       directory with all specified mpi_ranks.
# @param mpi_ranks:    A list of positive integers or "none"
# @param config_name:  A string describing the configuration to add
#                      this test to. Setting config_name to "" will
#                      add the test to the default configuration.
# @param work_dir:     A optional working directory argument where
#                      the test will be run.
#
# cmake-format: on
function(sparta_add_all_tests_to_config mpi_ranks config_name)
  set(__total_params 3)
  if(${ARGC} EQUAL ${__total_params})
    math(EXPR ARGC "${ARGC}-1")
    set(__work_dir "${ARGV${ARGC}}")
  else()
    set(__work_dir ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  # message("sparta_add_all_tests_to_config: __work_dir=${__work_dir}")

  file(
    GLOB __in_file_list
    LIST_DIRECTORIES false
    RELATIVE "${__work_dir}"
    CONFIGURE_DEPENDS ${__work_dir}/in.*)

  if(SPARTA_ENABLE_TESTING)
    if(SPARTA_DSMC_TESTING_PATH AND ${__work_dir} MATCHES
                                    "${SPARTA_DSMC_TESTING_ABSOLUTE_PATH}")
      # Skip disable list for tests in SPARTA_DSMC_TESTING_PATH
    else()
      list(REMOVE_ITEM __in_file_list ${SPARTA_DISABLED_TESTS})
    endif()
  endif()

  # message("__in_file_list=${__in_file_list}")
  if(BUILD_MPI)
    sparta_add_tests_to_config("${__in_file_list}" "${mpi_ranks}"
                               "${config_name}" ${__work_dir})
  else()
    sparta_add_tests_to_config("${__in_file_list}" "1" "${config_name}"
                               ${__work_dir})
  endif()
endfunction()

# Wrapper to sparta_add_all_tests_to_config
function(sparta_add_all_tests mpi_ranks)
  sparta_add_all_tests_to_config("${mpi_ranks}" "")
endfunction()
