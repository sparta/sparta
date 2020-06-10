if(PKG_MPI_STUBS)
  # ######### BEGIN TARGET_SPARTA_BUILD_MPI_STUBS ##########
  set(CMAKE_CXX_COMPILER "g++")
  set(CMAKE_CXX_FLAGS "-O -fPIC")
  set(CMAKE_AR "ar")

  # message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}:
  # CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

  set_source_files_properties(mpi.c PROPERTIES LANGUAGE CXX)
  add_library(${TARGET_SPARTA_BUILD_MPI} STATIC mpi.c)
  set_target_properties(${TARGET_SPARTA_BUILD_MPI} PROPERTIES LINKER_LANGUAGE
                                                              CXX)
  target_compile_options(${TARGET_SPARTA_BUILD_MPI}
                         PUBLIC ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})

  # Make include public so that targets which link against this can find the
  # includes
  target_include_directories(${TARGET_SPARTA_BUILD_MPI}
                             PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  install(TARGETS ${TARGET_SPARTA_BUILD_MPI}
          DESTINATION lib/sparta/${SPARTA_MACHINE}/pkgs)
  # ######### END   TARGET_SPARTA_BUILD_MPI_STUBS ##########
endif()
