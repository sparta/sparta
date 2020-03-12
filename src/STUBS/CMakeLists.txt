if(SPARTA_ENABLE_MPI_STUBS)
    ########## BEGIN TARGET_SPARTA_MPI_STUBS ##########
    add_library(${TARGET_SPARTA_MPI_STUBS} STATIC mpi.c)

    # Make include public so that targets which link against this can find the includes
    target_include_directories(${TARGET_SPARTA_MPI_STUBS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    ########## END   TARGET_SPARTA_MPI_STUBS ##########
endif()
