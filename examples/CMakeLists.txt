if(SPARTA_ENABLE_TESTING)
    include(${SPARTA_CMAKE_COMMON_DIR}/test/sparta_test_utils.cmake)

    file(GLOB SUITE_DIRS CONFIGURE_DEPENDS ${SPARTA_ENABLED_TEST_SUITES})
    file(COPY ${SUITE_DIRS} DESTINATION ${SPARTA_BINARY_DIR}/examples)

    foreach(SUITE ${SPARTA_ENABLED_TEST_SUITES})
        if(SPARTA_CTEST_CONFIGS)
            foreach(config ${SPARTA_CTEST_CONFIGS})
                # message("sparta_add_all_tests_to_config(${__DEFAULT_MPI_RANKS}
                # ${config})")
                sparta_add_all_tests_to_config("${__DEFAULT_MPI_RANKS}" "${config}" ${CMAKE_CURRENT_SOURCE_DIR}/${SUITE})
            endforeach()
        endif()
        sparta_add_all_tests_to_config("${__DEFAULT_MPI_RANKS}" "" ${CMAKE_CURRENT_SOURCE_DIR}/${SUITE})
    endforeach()
endif()
