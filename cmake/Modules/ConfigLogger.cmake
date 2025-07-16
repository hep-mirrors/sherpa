# Load CMake cache
get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(CACHE_VAR ${CACHE_VARS})
  get_property(CACHE_VAR_HELP CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
  # Capture compiler specifications
  if(CACHE_VAR_HELP STREQUAL "CXX compiler" OR
     CACHE_VAR_HELP STREQUAL "C compiler" OR
     CACHE_VAR_HELP STREQUAL "Fortran compiler")
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    if(CACHE_VAR_TYPE STREQUAL "STRING")
      set(CMDL_ARGS "${CMDL_ARGS} -D${CACHE_VAR}=\"${${CACHE_VAR}}\"")
    endif()
  endif()
  # Capture "enable"-type options with non-default values
  if(CACHE_VAR MATCHES "^.*_ENABLE_.*$")
    # First check the few disable-type flags
    if(CACHE_VAR MATCHES "^.*_ENABLE_EXAMPLES$" OR
       CACHE_VAR MATCHES "^.*_ENABLE_INTERNAL_PDFS$" OR
       CACHE_VAR MATCHES "^.*_ENABLE_LHAPDF$")
       if(${${CACHE_VAR}} STREQUAL "OFF")
         set(CMDL_ARGS "${CMDL_ARGS} -D${CACHE_VAR}=${${CACHE_VAR}}")
       endif()
    else() # then the enable-type flags
       if(${${CACHE_VAR}} STREQUAL "ON")
         set(CMDL_ARGS "${CMDL_ARGS} -D${CACHE_VAR}=${${CACHE_VAR}}")
       endif()
    endif()
  endif()
  # Capture path specifications
  if(CACHE_VAR_HELP STREQUAL "No help, variable specified on the command line." OR
     CACHE_VAR_HELP MATCHES "^.*for HepMC3.$")
    #get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    #if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
    #  set(CACHE_VAR_TYPE)
    #else()
    #  set(CACHE_VAR_TYPE ${CACHE_VAR_TYPE})
    #endif()
    #set(CMDL_ARGS "${CMDL_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
    set(CMDL_ARGS "${CMDL_ARGS} -D${CACHE_VAR}=\"${${CACHE_VAR}}\"")
  endif()
endforeach()
# Reconstruct the command line and write to config.log
file(WRITE "${CMAKE_BINARY_DIR}/config.log" "cmake${CMDL_ARGS} ${CMAKE_SOURCE_DIR}\n")
