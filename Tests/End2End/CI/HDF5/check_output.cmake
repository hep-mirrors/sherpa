if(NOT EXISTS "${OUTPUT_FILE}")
  message(FATAL_ERROR "Expected HDF5 output file was not created: ${OUTPUT_FILE}")
endif()

file(SIZE "${OUTPUT_FILE}" OUTPUT_SIZE)
if(OUTPUT_SIZE EQUAL 0)
  message(FATAL_ERROR "HDF5 output file is empty: ${OUTPUT_FILE}")
endif()