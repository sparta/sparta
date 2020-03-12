########## BEGIN TARGET_SPARTA ##########

########## BEGIN RECURSE COMMANDS ##########
# Required dependencies
# TODO: cmake,ccc,
add_subdirectory(STUBS)

# Optional Dependencies
## ADD_SUBDIRECTORY(FFT)
## ADD_SUBDIRECTORY(KOKKOS)
########## END   RECURSE COMMANDS ##########

########## START SETUP COMMANDS ##########
# GLOB together all source files
file(GLOB SPA_SRC_FILES
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        *.cpp
    )

# List out include directories
set(SPA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
########## END   SETUP COMMANDS ##########

########## START GENERATION COMMANDS ##########
add_executable(${TARGET_SPARTA} ${SPA_SRC_FILES})

target_include_directories(${TARGET_SPARTA} PRIVATE ${SPA_INCLUDE_DIRS})

target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_MPI_STUBS})

add_custom_command(TARGET ${TARGET_SPARTA} POST_BUILD COMMAND size ${TARGET_SPARTA})
########## END   GENERATION COMMANDS ##########

########## END   TARGET_SPARTA ##########
