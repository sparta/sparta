if(SPARTA_ENABLE_TESTING)
    set(SPARTA_ENABLED_TESTS 
            "ablation"
            "adapt"
        )
    file(GLOB test_files CONFIGURE_DEPENDS ${SPARTA_ENABLED_TESTS})
    file(COPY ${test_files} DESTINATION ${SPARTA_BINARY_DIR}/examples)

    foreach(TEST ${SPARTA_ENABLED_TESTS})
        message("Adding test: ${TEST}")
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${TEST})
    endforeach()
endif()
