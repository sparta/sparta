if(SPARTA_ENABLE_TESTING)
    include(${SPARTA_CMAKE_COMMON_DIR}/test/sparta_test_utils.cmake)

    set(SPARTA_ENABLED_TEST_SUITES
            "ablation"
            "adapt"
        )
    file(GLOB test_files CONFIGURE_DEPENDS ${SPARTA_ENABLED_TEST_SUITES})
    file(COPY ${test_files} DESTINATION ${SPARTA_BINARY_DIR}/examples)

    foreach(TEST ${SPARTA_ENABLED_TEST_SUITES})
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${TEST})
    endforeach()

    if(NOT SPARTA_DSMC_TESTING_PATH)
        message("Test correctness verfication is limited without SPARTA_DSMC_TESTING_PATH set!")
    endif()
endif()
