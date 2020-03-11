########## BEGIN RECURSE COMMANDS ##########
# Required dependencies
## add_subdirectory(STUBS)

# Optional Dependencies
## add_subdirectory(FFT)
## add_subdirectory(KOKKOS)
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
set(sparta_target spa_${__CMAKE_SPARTA_CONFIG})
########## END   SETUP COMMANDS ##########

########## START GENERATION COMMANDS ##########
add_executable(${sparta_target} ${SPA_SRC_FILES})
#add_executable(SPA ALIAS spa_${__CMAKE_SPARTA_CONFIG})
target_include_directories(${sparta_target} PRIVATE ${SPA_INCLUDE_DIRS})
########## END   GENERATION COMMANDS ##########
