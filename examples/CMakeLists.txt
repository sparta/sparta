if(SPARTA_ENABLE_TESTING)
    include(${SPARTA_CMAKE_COMMON_DIR}/test/sparta_test_utils.cmake)

    file(GLOB SUITE_DIRS CONFIGURE_DEPENDS ${SPARTA_ENABLED_TEST_SUITES})
    file(COPY ${SUITE_DIRS} DESTINATION ${SPARTA_BINARY_DIR}/examples)

    foreach(SUITE ${SPARTA_ENABLED_TEST_SUITES})
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${SUITE})
    endforeach()
endif()
