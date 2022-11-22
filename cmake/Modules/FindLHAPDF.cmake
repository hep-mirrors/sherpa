# - Try to find LHAPDF
# Defines:
#
#  LHAPDF_FOUND
#  LHAPDF_INCLUDE_DIR
#  LHAPDF_INCLUDE_DIRS (not cached)
#  LHAPDF_LIBRARY
#  LHAPDF_LIBRARIES (not cached)
#  LHAPDF_LIBRARY_DIR (not cached)

if (LHAPDF_ROOT_DIR OR LHAPDF_DIR OR (DEFINED ENV{LHAPDF_ROOT_DIR}) OR (DEFINED ENV{LHAPDF_DIR}) )
  set(LHAPDF_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (LHAPDF_ROOT_DIR)
    list (APPEND LHAPDF_SEARCH_DIRS "${LHAPDF_ROOT_DIR}" )
  endif()
  if (LHAPDF_DIR)
    list (APPEND LHAPDF_SEARCH_DIRS "${LHAPDF_DIR}" )
  endif()
  if (DEFINED EVN{LHAPDF_ROOT_DIR})
    list (APPEND LHAPDF_SEARCH_DIRS "$ENV{LHAPDF_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{LHAPDF_DIR})
    list (APPEND LHAPDF_SEARCH_DIRS "ENV{LHAPDF_DIR}" )
  endif()
endif()
if (LHAPDF_SEARCH_DIRS)
  find_program(LHAPDF_CONFIG_EXE NAMES lhapdf-config PATHS ${LHAPDF_SEARCH_DIRS} PATH_SUFFIXES bin NO_DEFAULT_PATH )
  find_path(LHAPDF_INCLUDE_DIR LHAPDF/LHAPDF.h PATHS ${LHAPDF_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(LHAPDF_LIBRARY NAMES LHAPDF PATHS ${LHAPDF_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_program(LHAPDF_CONFIG_EXE NAMES lhapdf-config PATHS ${LHAPDF_SEARCH_DIRS} PATH_SUFFIXES bin)
  find_path(LHAPDF_INCLUDE_DIR LHAPDF/LHAPDF.h PATH_SUFFIXES include)
  find_library(LHAPDF_LIBRARY NAMES LHAPDF PATHS_SUFFIXES lib lib64)
endif()
set(LHAPDF_VERSION 0.0.0)
if (LHAPDF_INCLUDE_DIR)
  if (EXISTS ${LHAPDF_INCLUDE_DIR}/LHAPDF/Version.h)
    file(STRINGS ${LHAPDF_INCLUDE_DIR}/LHAPDF/Version.h LHAPDF_VERSION_STRING_CONTENT REGEX "^#define[ ]+LHAPDF_VERSION[ ]+\"" )
    if (LHAPDF_VERSION_STRING_CONTENT)
      string(REGEX MATCH "[1234567890\.]+[a-zA-Z]*" LHAPDF_VERSION ${LHAPDF_VERSION_STRING_CONTENT})
    endif()
  endif()
endif()

set(LHAPDF_CONFIG_CPPFLAGS_STRING)
if (LHAPDF_CONFIG_EXE)
  execute_process(COMMAND ${LHAPDF_CONFIG_EXE} --cflags
                  OUTPUT_VARIABLE LHAPDF_CONFIG_CPPFLAGS_STRING
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

get_filename_component(LHAPDF_PATH ${LHAPDF_INCLUDE_DIR} DIRECTORY)
mark_as_advanced(LHAPDF_INCLUDE_DIR LHAPDF_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set LHAPDF_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LHAPDF DEFAULT_MSG LHAPDF_INCLUDE_DIR LHAPDF_LIBRARY LHAPDF_PATH LHAPDF_CONFIG_CPPFLAGS_STRING)


set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARY})
get_filename_component(LHAPDF_LIBRARY_DIR ${LHAPDF_LIBRARY} PATH)

set(LHAPDF_INCLUDE_DIRS ${LHAPDF_INCLUDE_DIR})



if(LHAPDF_FOUND AND NOT TARGET LHAPDF::LHAPDF)
    add_library(LHAPDF::LHAPDF UNKNOWN IMPORTED)
    set_target_properties(LHAPDF::LHAPDF PROPERTIES
        IMPORTED_LOCATION "${LHAPDF_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${LHAPDF_INCLUDE_DIRS}"
    )
endif()



mark_as_advanced(LHAPDF_FOUND)
