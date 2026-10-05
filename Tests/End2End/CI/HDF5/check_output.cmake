if(NOT EXISTS "${OUTPUT_FILE}")
  message(FATAL_ERROR "Expected HDF5 output file was not created: ${OUTPUT_FILE}")
endif()

file(SIZE "${OUTPUT_FILE}" OUTPUT_SIZE)
if(OUTPUT_SIZE EQUAL 0)
  message(FATAL_ERROR "HDF5 output file is empty: ${OUTPUT_FILE}")
endif()

# The content can only be checked if h5dump is available.
if(NOT H5DUMP)
  message(STATUS "h5dump not found, skip checking the HDF5 output content")
  return()
endif()

function(h5dump_output result)
  execute_process(
    COMMAND ${H5DUMP} ${ARGN} "${OUTPUT_FILE}"
    OUTPUT_VARIABLE output
    RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "h5dump ${ARGN} failed on ${OUTPUT_FILE}")
  endif()
  set(${result} "${output}" PARENT_SCOPE)
endfunction()

function(expect output regex description)
  if(NOT output MATCHES "${regex}")
    message(FATAL_ERROR "HDF5 output: ${description}\n${output}")
  endif()
endfunction()

h5dump_output(version -d version)
expect("${version}" "\\(0\\): 2, 7, 0" "unexpected LHEH5 format version")

h5dump_output(attributes -A)
foreach(column EXTRA__MEWeight EXTRA__PSWeight numRandomNumbers)
  expect("${attributes}" "\"${column}\"" "column ${column} missing")
endforeach()
foreach(attribute generatorName generatorVersion phasespaceGenerator
                  description)
  expect("${attributes}" "ATTRIBUTE \"${attribute}\""
         "init attribute ${attribute} missing")
endforeach()
expect("${attributes}" "DATASET \"generatedResult\""
       "generatedResult dataset missing")
if(attributes MATCHES "ATTRIBUTE \"(events|ctevents)\"")
  message(FATAL_ERROR "HDF5 output: column names not stored as \"properties\"")
endif()

# aqed of the first event, which must be set
h5dump_output(aqed -d events -s 0,7 -c 1,1)
expect("${aqed}" "\\(0,7\\): 0\\.00" "alpha_QED not written")
