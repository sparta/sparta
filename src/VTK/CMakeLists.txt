if(PKG_VTK)
  # ######### BEGIN TARGET_SPARTA_PKG_VTK ##########

  file(
    GLOB SPARTA_PKG_VTK_SRC_FILES
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS *.cpp *.hpp *.h)

  include(${SPARTA_CMAKE_COMMON_DIR}/set/style_file_glob.cmake)

  add_library(${TARGET_SPARTA_PKG_VTK} STATIC ${SPARTA_PKG_VTK_SRC_FILES})
  set_target_properties(${TARGET_SPARTA_PKG_VTK} PROPERTIES LINKER_LANGUAGE CXX)

  # Get MPI header files
  target_link_libraries(${TARGET_SPARTA_PKG_VTK} LINK_PRIVATE
                        ${TARGET_SPARTA_BUILD_MPI})

  # Get Kokkos header files
  target_link_libraries(${TARGET_SPARTA_PKG_VTK} LINK_PRIVATE
                        ${TARGET_SPARTA_BUILD_KOKKOS})

  # Get VTK header files and libraries
  target_link_libraries(${TARGET_SPARTA_PKG_VTK} LINK_PUBLIC
                        ${TARGET_SPARTA_BUILD_VTK})

  # Make include public so that targets which link against this can find the
  # includes
  target_include_directories(${TARGET_SPARTA_PKG_VTK}
                             PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  # Add include dependencies for building TARGET_SPARTA_PKG_VTK
  target_include_directories(${TARGET_SPARTA_PKG_VTK}
                             PRIVATE ${SPA_INCLUDE_DIRS})

  # VTK 9 module system requires autoinit so the IO writer object factories
  # register at load time; without this the XML/legacy writers produce no output
  vtk_module_autoinit(TARGETS ${TARGET_SPARTA_PKG_VTK}
                      MODULES ${TARGET_SPARTA_BUILD_VTK})

  install(TARGETS ${TARGET_SPARTA_PKG_VTK}
          DESTINATION lib/sparta/${SPARTA_MACHINE}/pkgs)
  # ######### END   TARGET_SPARTA_PKG_VTK ##########
endif()
