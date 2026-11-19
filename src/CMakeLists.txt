# ######### BEGIN TARGET_SPARTA ##########

# ######### START SETUP COMMANDS ##########
# GLOB together all source files
file(
  GLOB SPA_SRC_FILES
  LIST_DIRECTORIES false
  CONFIGURE_DEPENDS *.cpp *.hpp *.h)
list(REMOVE_ITEM SPA_SRC_FILES "main.cpp")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# List out include directories
set(SPA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_BINARY_DIR}/include)
set(SPARTA_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
# ######### END   SETUP COMMANDS ##########

# ######### BEGIN PROCESSING PACKAGES - RECURSE COMMANDS ##########
add_subdirectory(${SPARTA_SRC_DIR}/KOKKOS)
add_subdirectory(${SPARTA_SRC_DIR}/FFT)
add_subdirectory(${SPARTA_SRC_DIR}/STUBS)
add_subdirectory(${SPARTA_SRC_DIR}/PYTHON)
add_subdirectory(${SPARTA_SRC_DIR}/VTK)
# message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}:
# CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
# ######### END   PROCESSING PACKAGES - RECURSE COMMANDS ##########

# ######### START GENERATION COMMANDS ##########
# message(VERBOSE "TARGET_SPARTA_BUILD_TPLS: ${TARGET_SPARTA_BUILD_TPLS}")
# message(VERBOSE "TARGET_SPARTA_PKGS: ${TARGET_SPARTA_PKGS}") message(VERBOSE
# "SPA_INCLUDE_DIRS: ${SPA_INCLUDE_DIRS}")
add_library(${TARGET_SPARTA_LIB} ${SPA_SRC_FILES})
set_target_properties(
  ${TARGET_SPARTA_LIB} PROPERTIES VERSION ${SPARTA_BUILD_VER} SOVERSION
                                                              ${SPARTA_API_VER})
target_include_directories(${TARGET_SPARTA_LIB} PRIVATE ${SPA_INCLUDE_DIRS})
target_link_libraries(${TARGET_SPARTA_LIB} LINK_PUBLIC
                      ${TARGET_SPARTA_BUILD_TPLS})
target_link_libraries(${TARGET_SPARTA_LIB} LINK_PUBLIC ${TARGET_SPARTA_PKGS})
sparta_kokkos_host_only_sources(${TARGET_SPARTA_LIB})

# Generate spartagitversion.h with the git commit and branch at every build,
# as LAMMPS does.  Only sparta.cpp includes it, and it is only replaced when
# its contents change.
find_package(Git QUIET)
add_custom_target(
  ${TARGET_SPARTA_LIB}_GITVERSION
  COMMAND
    ${CMAKE_COMMAND} -DSPARTA_DIR="${SPARTA_SRC_DIR}/.."
    -DSPARTA_GIT_HEADER="${SPARTA_BINARY_DIR}/include/spartagitversion.h"
    -DGIT_FOUND="${GIT_FOUND}" -DGIT_EXECUTABLE="${GIT_EXECUTABLE}" -P
    ${SPARTA_CMAKE_DIR}/Modules/generate_spartagitversion.cmake
  BYPRODUCTS ${SPARTA_BINARY_DIR}/include/spartagitversion.h)
add_dependencies(${TARGET_SPARTA_LIB} ${TARGET_SPARTA_LIB}_GITVERSION)
set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sparta.cpp APPEND
             PROPERTY COMPILE_DEFINITIONS SPARTA_GIT_VERSION_HEADER)

include(${SPARTA_CMAKE_COMMON_DIR}/set/style_file_glob.cmake)

add_custom_target(
  ${TARGET_SPARTA_LIB}_PRE_BUILD
  COMMAND sh ${SPARTA_CMAKE_DIR}/make_style.sh ${SPARTA_BINARY_DIR}/include
          ${CMAKE_CURRENT_SOURCE_DIR})

add_dependencies(${TARGET_SPARTA_LIB} ${TARGET_SPARTA_LIB}_PRE_BUILD)
add_custom_command(
  TARGET ${TARGET_SPARTA_LIB}
  POST_BUILD
  COMMAND "${CMAKE_COMMAND}" -E create_symlink ${TARGET_SPARTA_LIB_NAME}
          ${TARGET_SPARTA_LIB_LINK})
install(TARGETS ${TARGET_SPARTA_LIB} DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_SPARTA_LIB_LINK}
        DESTINATION lib)

add_executable(${TARGET_SPARTA} main.cpp)
target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_LIB})
sparta_kokkos_host_only_sources(${TARGET_SPARTA})
# print the binary size after linking (informational). Skip when cross-
# compiling: the command uses the bare target name, but on Windows the file is
# spa_.exe, so `size spa_` fails the build during a MinGW64 cross-compile.
if(NOT CMAKE_CROSSCOMPILING)
  add_custom_command(
    TARGET ${TARGET_SPARTA}
    POST_BUILD
    COMMAND size ${TARGET_SPARTA})
endif()
install(TARGETS ${TARGET_SPARTA} DESTINATION bin)
# ######### END   GENERATION COMMANDS ##########

# ######### END   TARGET_SPARTA ##########
