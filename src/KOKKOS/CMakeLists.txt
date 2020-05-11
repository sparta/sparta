if(PKG_KOKKOS)
    ########## BEGIN TARGET_SPARTA_PKG_KOKKOS ##########
    #set(CMAKE_CXX_COMPILER "c++")
    message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}: CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
    
    if(BUILD_KOKKOS_TPL)
        add_subdirectory(${SPARTA_TPL_DIR}/kokkos ${CMAKE_CURRENT_BINARY_DIR}/lib/kokkos)
        
        message(VERBOSE "Kokkos_INCLUDE_DIRS_RET: ${Kokkos_INCLUDE_DIRS_RET}")
        message(VERBOSE "TARGET_SPARTA_BUILD_KOKKOS_TPL: ${TARGET_SPARTA_BUILD_KOKKOS_TPL}")
        get_target_property(TARGET_SPARTA_BUILD_KOKKOS_TPL_INCLUDE_DIRECTORUES ${TARGET_SPARTA_BUILD_KOKKOS_TPL} INTERFACE_INCLUDE_DIRECTORIES)
        message(VERBOSE "TARGET_SPARTA_BUILD_KOKKOS_TPL_INCLUDE_DIRECTORUES: ${TARGET_SPARTA_BUILD_KOKKOS_TPL_INCLUDE_DIRECTORUES}")
    endif()

    file(GLOB SPARTA_PKG_KOKKOS_SRC_FILES
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        *.cpp *.hpp *.h
    )

    add_library(${TARGET_SPARTA_PKG_KOKKOS} STATIC ${SPARTA_PKG_KOKKOS_SRC_FILES})
    set_target_properties(${TARGET_SPARTA_PKG_KOKKOS} PROPERTIES LINKER_LANGUAGE CXX)

    target_link_libraries(${TARGET_SPARTA_PKG_KOKKOS} LINK_PRIVATE ${TARGET_SPARTA_BUILD_MPI_TPL})
    
    # This should handle resolving include directories, but these include paths
    # are actually returned via Kokkos_INCLUDE_DIRS_RET.
    target_link_libraries(${TARGET_SPARTA_PKG_KOKKOS} LINK_PRIVATE ${TARGET_SPARTA_BUILD_KOKKOS_TPL})
    
    set_target_properties(${TARGET_SPARTA_PKG_KOKKOS} PROPERTIES COMPILE_FLAGS "${SPARTA_CXX_COMPILE_FLAGS} ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS}")
    
    # Make include public so that targets which links against this can find the includes
    target_include_directories(${TARGET_SPARTA_PKG_KOKKOS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${TARGET_SPARTA_PKG_KOKKOS} PUBLIC ${Kokkos_INCLUDE_DIRS_RET})
    
    # Add include dependencies for building TARGET_SPARTA_PKG_KOKKOS
    target_include_directories(${TARGET_SPARTA_PKG_KOKKOS} PRIVATE ${SPA_INCLUDE_DIRS})
    ########## END   TARGET_SPARTA_PKG_KOKKOS ##########
endif()
