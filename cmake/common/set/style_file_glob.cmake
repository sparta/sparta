file(GLOB style_files "*.h")
foreach(file ${style_files})
  configure_file(${file} ${SPARTA_BINARY_DIR}/include COPYONLY)
endforeach()
