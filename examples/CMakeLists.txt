if(SPARTA_ENABLE_TESTING)
    include(${SPARTA_CMAKE_COMMON_DIR}/test/sparta_test_utils.cmake)

    set(SPARTA_ENABLED_TEST_SUITES
            "ablation"
            "adapt"
        )

    list(APPEND __DEFAULT_MPI_RANKS "1")
    list(APPEND __DEFAULT_MPI_RANKS "4")

    file(GLOB test_files CONFIGURE_DEPENDS ${SPARTA_ENABLED_TEST_SUITES})
    file(COPY ${test_files} DESTINATION ${SPARTA_BINARY_DIR}/examples)

    foreach(TEST ${SPARTA_ENABLED_TEST_SUITES})
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${TEST})
    endforeach()
endif()
