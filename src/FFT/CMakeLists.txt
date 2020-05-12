if(PKG_FFT)
    ########## BEGIN TARGET_SPARTA_PKG_FFT ##########
    message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}: CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

    file(GLOB SPARTA_PKG_FFT_SRC_FILES
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        *.cpp *.hpp *.h
    )

    add_library(${TARGET_SPARTA_PKG_FFT} STATIC ${SPARTA_PKG_FFT_SRC_FILES})
    set_target_properties(${TARGET_SPARTA_PKG_FFT} PROPERTIES LINKER_LANGUAGE CXX)

    # Get MPI header files
    target_link_libraries(${TARGET_SPARTA_PKG_FFT} LINK_PRIVATE ${TARGET_SPARTA_BUILD_MPI_TPL})
    
    # Make include public so that targets which links against this can find the includes
    target_include_directories(${TARGET_SPARTA_PKG_FFT} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    
    # Add include dependencies for building TARGET_SPARTA_PKG_FFT
    target_include_directories(${TARGET_SPARTA_PKG_FFT} PRIVATE ${SPA_INCLUDE_DIRS})
    ########## END   TARGET_SPARTA_PKG_FFT ##########
endif()
