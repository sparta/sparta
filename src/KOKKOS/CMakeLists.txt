if(PKG_KOKKOS)
  # ######### BEGIN TARGET_SPARTA_PKG_KOKKOS ##########
  # message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}:
  # CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

  if(BUILD_KOKKOS)
    include(${SPARTA_CMAKE_COMMON_DIR}/set/kokkos_cmake_defaults.cmake)

    # message(VERBOSE "Checking Kokkos_ENABLE_CUDA: ${Kokkos_ENABLE_CUDA}")
    if(Kokkos_ENABLE_CUDA)
      include(${SPARTA_CMAKE_COMMON_DIR}/set/kokkos_cuda_cmake_defaults.cmake)
      set(ENV{MPICH_CXX} ${SPARTA_TPL_DIR}/kokkos/bin/nvcc_wrapper)
      set(ENV{OMPI_CXX} $ENV{MPICH_CXX})
      # message(VERBOSE "MPICH_CXX: $ENV{MPICH_CXX}") message(VERBOSE "OMPI_CXX:
      # $ENV{OMPI_CXX}")
    endif()

    add_subdirectory(${SPARTA_TPL_DIR}/kokkos
                     ${CMAKE_CURRENT_BINARY_DIR}/lib/kokkos)

    # message(VERBOSE "Kokkos_INCLUDE_DIRS_RET: ${Kokkos_INCLUDE_DIRS_RET}")
    # message(VERBOSE "TARGET_SPARTA_BUILD_KOKKOS:
    # ${TARGET_SPARTA_BUILD_KOKKOS}")
    get_target_property(
      TARGET_SPARTA_BUILD_KOKKOS_INCLUDE_DIRECTORUES
      ${TARGET_SPARTA_BUILD_KOKKOS} INTERFACE_INCLUDE_DIRECTORIES)
    # message(VERBOSE "TARGET_SPARTA_BUILD_KOKKOS_INCLUDE_DIRECTORUES:
    # ${TARGET_SPARTA_BUILD_KOKKOS_INCLUDE_DIRECTORUES}")
  endif()

  file(
    GLOB SPARTA_PKG_KOKKOS_SRC_FILES
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS *.cpp *.hpp *.h)

  include(${SPARTA_CMAKE_COMMON_DIR}/set/style_file_glob.cmake)

  add_library(${TARGET_SPARTA_PKG_KOKKOS} STATIC ${SPARTA_PKG_KOKKOS_SRC_FILES})
  set_target_properties(${TARGET_SPARTA_PKG_KOKKOS} PROPERTIES LINKER_LANGUAGE
                                                               CXX)

  target_link_libraries(${TARGET_SPARTA_PKG_KOKKOS} LINK_PRIVATE
                        ${TARGET_SPARTA_BUILD_MPI})

  # This should handle resolving include directories when kokkos is updated.
  target_link_libraries(${TARGET_SPARTA_PKG_KOKKOS} LINK_PRIVATE
                        ${TARGET_SPARTA_BUILD_KOKKOS})

  # Make include public so that targets which links against this can find the
  # includes
  target_include_directories(${TARGET_SPARTA_PKG_KOKKOS}
                             PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${TARGET_SPARTA_PKG_KOKKOS} SYSTEM
                             PUBLIC ${Kokkos_INCLUDE_DIRS_RET})

  # Add include dependencies for building TARGET_SPARTA_PKG_KOKKOS
  target_include_directories(${TARGET_SPARTA_PKG_KOKKOS}
                             PRIVATE ${SPA_INCLUDE_DIRS})

  install(TARGETS ${TARGET_SPARTA_PKG_KOKKOS}
          DESTINATION lib/sparta/${SPARTA_MACHINE}/pkgs)
  # ######### END   TARGET_SPARTA_PKG_KOKKOS ##########
endif()
