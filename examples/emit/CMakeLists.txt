if(SPARTA_CTEST_CONFIGS)
    foreach(config ${SPARTA_CTEST_CONFIGS})
        # message("sparta_add_all_tests_to_config(${__DEFAULT_MPI_RANKS} ${config})")
        sparta_add_all_tests_to_config("${__DEFAULT_MPI_RANKS}" "${config}")
    endforeach()
endif()

sparta_add_all_tests_to_config("${__DEFAULT_MPI_RANKS}" "")
