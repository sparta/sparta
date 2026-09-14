if(PKG_KOKKOS)
  # ######### BEGIN TARGET_SPARTA_PKG_KOKKOS ##########
  # message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}:
  # CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

  ########################################################################
  # consistency checks and Kokkos options/settings required by SPARTA

  # Kokkos 5.2 still provides the legacy view implementation, which computes
  # view offsets in 64-bit arithmetic; the newer one does not, and can index a
  # view of more than 2^31 entries incorrectly.  Keep the legacy views until the
  # KOKKOS package is ready for the newer ones (Kokkos 5.3 removes them), and do
  # not let the build be talked out of it.  The legacy implementation is
  # deprecated in 5.2, so silence the warnings about it.
  set(Kokkos_ENABLE_IMPL_VIEW_LEGACY ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_DEPRECATION_WARNINGS OFF CACHE BOOL "" FORCE)
  mark_as_advanced(Kokkos_ENABLE_IMPL_VIEW_LEGACY Kokkos_ENABLE_DEPRECATION_WARNINGS)

  # Kokkos removed Kokkos_ENABLE_CUDA_UVM in 5.2 (deprecated since 4.0); CMake
  # would otherwise accept the stale setting and silently do nothing with it
  if(Kokkos_ENABLE_CUDA_UVM)
    message(FATAL_ERROR "The option Kokkos_ENABLE_CUDA_UVM is no longer supported by Kokkos. "
      "Use -D Kokkos_ENABLE_IMPL_CUDA_UNIFIED_MEMORY=on instead, which requires CUDA 12.2 or later")
  endif()

  if(Kokkos_ENABLE_CUDA)
    option(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC "CUDA asynchronous malloc support" OFF)
    mark_as_advanced(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC)
    if(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC)
      message(STATUS "KOKKOS: CUDA malloc async support enabled")
    else()
      message(STATUS "KOKKOS: CUDA malloc async support disabled")
    endif()

    option(Kokkos_ENABLE_IMPL_CUDA_UNIFIED_MEMORY "Enable unified memory with CUDA" OFF)
    mark_as_advanced(Kokkos_ENABLE_IMPL_CUDA_UNIFIED_MEMORY)
    if(Kokkos_ENABLE_IMPL_CUDA_UNIFIED_MEMORY)
      message(STATUS "KOKKOS: CUDA unified memory enabled")
    else()
      message(STATUS "KOKKOS: CUDA unified memory disabled")
    endif()
  endif()

  if(Kokkos_ENABLE_HIP)
    option(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS "Enable multiple kernel instantiations with HIP" ON)
    mark_as_advanced(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS)
    option(Kokkos_ENABLE_ROCTHRUST "Use RoCThrust library" ON)
    mark_as_advanced(Kokkos_ENABLE_ROCTHRUST)

    if(Kokkos_ARCH_AMD_GFX942 OR Kokkos_ARCH_AMD_GFX940)
      option(Kokkos_ENABLE_IMPL_HIP_UNIFIED_MEMORY "Enable unified memory with HIP" ON)
      mark_as_advanced(Kokkos_ENABLE_IMPL_HIP_UNIFIED_MEMORY)
    endif()
  endif()

  if(Kokkos_ENABLE_SERIAL)
    if(NOT (Kokkos_ENABLE_OPENMP OR Kokkos_ENABLE_THREADS OR
      Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL
      OR Kokkos_ENABLE_OPENMPTARGET))
    option(Kokkos_ENABLE_ATOMICS_BYPASS "Disable atomics for Kokkos Serial Backend" ON)
    mark_as_advanced(Kokkos_ENABLE_ATOMICS_BYPASS)
    endif()
  endif()
  ########################################################################

  if(BUILD_KOKKOS)
    if(USE_EXTERNAL_KOKKOS)
      find_package(Kokkos REQUIRED CONFIG)
      set(TARGET_SPARTA_BUILD_KOKKOS Kokkos::kokkos)
      get_target_property(Kokkos_INCLUDE_DIRS_RET ${TARGET_SPARTA_BUILD_KOKKOS} INTERFACE_INCLUDE_DIRECTORIES)
    else()
      include(${SPARTA_CMAKE_COMMON_DIR}/set/kokkos_cmake_defaults.cmake)

      # message(VERBOSE "Checking Kokkos_ENABLE_CUDA: ${Kokkos_ENABLE_CUDA}")
      if(Kokkos_ENABLE_CUDA)
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
  endif()

  file(
    GLOB SPARTA_PKG_KOKKOS_SRC_FILES
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS *.cpp *.hpp *.h)

  if (NOT PKG_FFT)
    list(FILTER SPARTA_PKG_KOKKOS_SRC_FILES EXCLUDE REGEX ".*fft.*kokkos.*")
    list(FILTER SPARTA_PKG_KOKKOS_SRC_FILES EXCLUDE REGEX ".*pack.*kokkos.*")
    list(FILTER SPARTA_PKG_KOKKOS_SRC_FILES EXCLUDE REGEX ".*remap.*kokkos.*")
    list(REMOVE_ITEM style_files kokkos_base_fft.h)
  endif()

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

  message(VERBOSE "Kokkos_INCLUDE_DIRS_RET: ${Kokkos_INCLUDE_DIRS_RET}")  
  message(VERBOSE "TARGET_SPARTA_PKG_KOKKOS: ${TARGET_SPARTA_PKG_KOKKOS}")  

  # Add include dependencies for building TARGET_SPARTA_PKG_KOKKOS
  target_include_directories(${TARGET_SPARTA_PKG_KOKKOS}
                             PRIVATE ${SPA_INCLUDE_DIRS})

  install(TARGETS ${TARGET_SPARTA_PKG_KOKKOS}
          DESTINATION lib/sparta/${SPARTA_MACHINE}/pkgs)
  # ######### END   TARGET_SPARTA_PKG_KOKKOS ##########
endif()
