########## BEGIN TARGET_SPARTA ##########

########## START SETUP COMMANDS ##########
# GLOB together all source files
file(GLOB SPA_SRC_FILES
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        *.cpp *.hpp *.h
    )

# List out include directories
set(SPA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
message(VERBOSE "SPA_INCLUDE_DIRS=${SPA_INCLUDE_DIRS}")
########## END   SETUP COMMANDS ##########

########## BEGIN PROCESSING PACKAGES - RECURSE COMMANDS ##########
add_subdirectory(${SPARTA_PKG_DIR}/FFT)
add_subdirectory(${SPARTA_PKG_DIR}/KOKKOS)
add_subdirectory(${SPARTA_PKG_DIR}/STUBS)
message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}: CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
########## END   PROCESSING PACKAGES - RECURSE COMMANDS ##########

########## START GENERATION COMMANDS ##########
add_executable(${TARGET_SPARTA} ${SPA_SRC_FILES})

target_include_directories(${TARGET_SPARTA} PRIVATE ${SPA_INCLUDE_DIRS})

target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_BUILD_TPLS})
target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_PKGS})

# TODO: Eliminate PRE_BUILD target
add_custom_target(${TARGET_SPARTA}_PRE_BUILD COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/Make.sh style)
add_dependencies(${TARGET_SPARTA} ${TARGET_SPARTA}_PRE_BUILD)

add_custom_command(TARGET ${TARGET_SPARTA} POST_BUILD COMMAND size ${TARGET_SPARTA})
########## END   GENERATION COMMANDS ##########

install(TARGETS ${TARGET_SPARTA} DESTINATION bin)
########## END   TARGET_SPARTA ##########
