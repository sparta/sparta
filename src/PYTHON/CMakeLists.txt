if(PKG_PYTHON)
  # ######### BEGIN TARGET_SPARTA_PKG_PYTHON ##########
  # message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}:
  # CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

  file(
    GLOB SPARTA_PKG_PYTHON_SRC_FILES
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS *.cpp *.hpp *.h)

  include(${SPARTA_CMAKE_COMMON_DIR}/set/style_file_glob.cmake)

  add_library(${TARGET_SPARTA_PKG_PYTHON} STATIC ${SPARTA_PKG_PYTHON_SRC_FILES})
  set_target_properties(${TARGET_SPARTA_PKG_PYTHON} PROPERTIES LINKER_LANGUAGE CXX)

  # Get MPI header files
  target_link_libraries(${TARGET_SPARTA_PKG_PYTHON} LINK_PRIVATE
                        ${TARGET_SPARTA_BUILD_MPI})

  # Get Kokkos header files
  target_link_libraries(${TARGET_SPARTA_PKG_PYTHON} LINK_PRIVATE
                        ${TARGET_SPARTA_BUILD_KOKKOS})

  # Get PYTHON header files
  if(PYTHON)
    target_link_libraries(${TARGET_SPARTA_PKG_PYTHON} LINK_PRIVATE
                          ${TARGET_SPARTA_BUILD_PYTHON})
  endif()

  # Make include public so that targets which links against this can find the
  # includes
  target_include_directories(${TARGET_SPARTA_PKG_PYTHON}
                             PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  # Add include dependencies for building TARGET_SPARTA_PKG_PYTHON
  target_include_directories(${TARGET_SPARTA_PKG_PYTHON}
                             PRIVATE ${SPA_INCLUDE_DIRS})

  install(TARGETS ${TARGET_SPARTA_PKG_PYTHON}
          DESTINATION lib/sparta/${SPARTA_MACHINE}/pkgs)
  # ######### END   TARGET_SPARTA_PKG_PYTHON ##########
endif()
