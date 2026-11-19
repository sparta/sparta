file(GLOB style_files "*.h")

if (PKG_KOKKOS)
  if (NOT PKG_FFT)
    list(FILTER style_files EXCLUDE REGEX "/[^/]*(fft|pack|remap)[^/]*kokkos[^/]*$")
    list(FILTER style_files EXCLUDE REGEX "/kokkos_base_fft\\.h$")
  endif()
endif()

foreach(file ${style_files})
  configure_file(${file} ${SPARTA_BINARY_DIR}/include COPYONLY)
endforeach()
