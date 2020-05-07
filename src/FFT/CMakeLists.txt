if(PACKAGE_FFT)
    ########## BEGIN TARGET_SPARTA_FFT ##########
    #set(CMAKE_CXX_COMPILER "c++")

    message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}: CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

    file(GLOB SPARTA_PACKAGE_FFT_SRC_FILES
        LIST_DIRECTORIES false
        CONFIGURE_DEPENDS
        *.cpp
    )

    add_library(${TARGET_SPARTA_PACKAGE_FFT} STATIC ${SPARTA_PACKAGE_FFT_SRC_FILES})
    set_target_properties(${TARGET_SPARTA_PACKAGE_FFT} PROPERTIES LINKER_LANGUAGE CXX)

    # Get MPI header files
    target_link_libraries(${TARGET_SPARTA_PACKAGE_FFT} LINK_PRIVATE ${TARGET_SPARTA_BUILD_MPI_TPL})
    
    set_target_properties(${TARGET_SPARTA_PACKAGE_FFT} PROPERTIES COMPILE_FLAGS "${SPARTA_CXX_COMPILE_FLAGS} ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS}")
    
    # Make include public so that targets which links against this can find the includes
    target_include_directories(${TARGET_SPARTA_PACKAGE_FFT} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${TARGET_SPARTA_PACKAGE_FFT} PRIVATE ${SPA_INCLUDE_DIRS})
    ########## END   TARGET_SPARTA_FFT ##########
endif()

#cd /Users/eharvey/dev/SPARTA.base/sparta/build/src && /usr/bin/c++   -I/Users/eharvey/dev/SPARTA.base/sparta/src -I/Users/eharvey/dev/SPARTA.base/sparta/build/src -I/Users/eharvey/dev/SPARTA.base/sparta/src/STUBS  -O -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk   -fPIC -DSPARTA_UNORDERED_MAP -DSPARTA_GZIP -o CMakeFiles/spa_.dir/collide.cpp.o -c /Users/eharvey/dev/SPARTA.base/sparta/src/collide.cpp
#/Users/eharvey/dev/SPARTA.base/sparta/build/src/FFT && /usr/bin/c++   -I/Users/eharvey/dev/SPARTA.base/sparta/src/FFT -I/Users/eharvey/dev/SPARTA.base/sparta/src -I/Users/eharvey/dev/SPARTA.base/sparta/build/src -I/Users/eharvey/dev/SPARTA.base/sparta/src/STUBS  -O -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk   -o CMakeFiles/fft.dir/compute_fft_grid.cpp.o -c /Users/eharvey/dev/SPARTA.base/sparta/src/FFT/compute_fft_grid.cpp
