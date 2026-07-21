if(PKG_MPI_STUBS)
  # ######### BEGIN TARGET_SPARTA_BUILD_MPI_STUBS ##########
  # Force the plain host compiler/archiver for the trivial MPI stubs when
  # building natively (e.g. so a Kokkos/CUDA wrapper compiler is not used for
  # them). Do NOT do this when cross-compiling: the stubs must be built with
  # the same cross toolchain as the rest of SPARTA, otherwise they produce
  # host-format objects that the target linker cannot use and the MPI_* symbols
  # go missing from libsparta (e.g. the MinGW64 Windows DLL link).
  if(NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_CXX_COMPILER "g++")
    set(CMAKE_CXX_FLAGS "-O -fPIC")
    set(CMAKE_AR "ar")
  endif()

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
