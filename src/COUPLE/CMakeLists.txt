if(PKG_COUPLE)
  # ######### BEGIN TARGET_SPARTA_PKG_COUPLE ##########
  # message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}:
  # CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

  file(
    GLOB SPARTA_PKG_COUPLE_SRC_FILES
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS *.cpp *.hpp *.h)

  include(${SPARTA_CMAKE_COMMON_DIR}/set/style_file_glob.cmake)

  add_library(${TARGET_SPARTA_PKG_COUPLE} STATIC ${SPARTA_PKG_COUPLE_SRC_FILES})
  set_target_properties(${TARGET_SPARTA_PKG_COUPLE} PROPERTIES LINKER_LANGUAGE CXX)

  # Get MPI header files
  target_link_libraries(${TARGET_SPARTA_PKG_COUPLE} LINK_PRIVATE
                        ${TARGET_SPARTA_BUILD_MPI})

  find_library(LAMMPS_LIBRARY lammps ${SPARTA_SRC_DIR}/../lib/lammps/build)
  if(NOT LAMMPS_LIBRARY)
    message(FATAL_ERROR "LAMMPS must be built in sparta/lib/lammps")
  endif()

  target_link_libraries(${TARGET_SPARTA_PKG_COUPLE} LINK_PRIVATE
                      ${LAMMPS_LIBRARY})

  # Make include public so that targets which links against this can find the
  # includes
  target_include_directories(${TARGET_SPARTA_PKG_COUPLE}
                             PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  # Add include dependencies for building TARGET_SPARTA_PKG_COUPLE
  target_include_directories(${TARGET_SPARTA_PKG_COUPLE}
                             PRIVATE ${SPA_INCLUDE_DIRS})

  target_include_directories(${TARGET_SPARTA_PKG_COUPLE}
                             PUBLIC ${SPARTA_SRC_DIR}/../lib/lammps/src)

  install(TARGETS ${TARGET_SPARTA_PKG_COUPLE}
          DESTINATION lib/sparta/${SPARTA_MACHINE}/pkgs)
  # ######### END   TARGET_SPARTA_PKG_COUPLE ##########
endif()
