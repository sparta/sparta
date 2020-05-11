########## BEGIN TARGET_SPARTA ##########

########## START SETUP COMMANDS ##########
# GLOB together all source files
file(GLOB SPA_SRC_FILES
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        *.cpp
    )

# List out include directories
set(SPA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
message(VERBOSE "SPA_INCLUDE_DIRS=${SPA_INCLUDE_DIRS}")
########## END   SETUP COMMANDS ##########

########## BEGIN RECURSE COMMANDS ##########
# Required dependencies
# TODO: cmake,ccc,
add_subdirectory(FFT)
add_subdirectory(KOKKOS)
add_subdirectory(STUBS)
message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}: CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

# Optional Dependencies
## ADD_SUBDIRECTORY(FFT)
## ADD_SUBDIRECTORY(KOKKOS)
########## END   RECURSE COMMANDS ##########

########## START GENERATION COMMANDS ##########
add_executable(${TARGET_SPARTA} ${SPA_SRC_FILES})

target_include_directories(${TARGET_SPARTA} PRIVATE ${SPA_INCLUDE_DIRS})

target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_BUILD_TPLS})
target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_PKGS})

add_custom_target(${TARGET_SPARTA}_PRE_BUILD COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/Make.sh style)

add_dependencies(${TARGET_SPARTA} ${TARGET_SPARTA}_PRE_BUILD)

add_custom_command(TARGET ${TARGET_SPARTA} POST_BUILD COMMAND size ${TARGET_SPARTA})
########## END   GENERATION COMMANDS ##########

########## END   TARGET_SPARTA ##########
