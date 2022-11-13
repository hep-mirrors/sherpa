# - Try to find GOSAM
# Defines:
#
#  GOSAM_FOUND
#  GOSAM_PREFIX
if (GOSAM_ROOT_DIR OR GOSAM_DIR OR (DEFINED ENV{GOSAM_ROOT_DIR}) OR (DEFINED ENV{GOSAM_DIR}) )
  set(GOSAM_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (GOSAM_ROOT_DIR)
    list (APPEND GOSAM_SEARCH_DIRS "${GOSAM_ROOT_DIR}" )
  endif()
  if (GOSAM_DIR)
    list (APPEND GOSAM_SEARCH_DIRS "${GOSAM_DIR}" )
  endif()
  if (DEFINED EVN{GOSAM_ROOT_DIR})
    list (APPEND GOSAM_SEARCH_DIRS "$ENV{GOSAM_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{GOSAM_DIR})
    list (APPEND GOSAM_SEARCH_DIRS "ENV{GOSAM_DIR}" )
  endif()
endif()
set(GOSAM_VERSION Unknown)
find_program(GOSAM NAMES gosam gosam.py )
if (GOSAM)
execute_process(COMMAND ${GOSAM} --version
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE GOSAM_VERSION_ALL 
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
string(SUBSTRING "${GOSAM_VERSION_ALL}" 6 5 GOSAM_VERSION)
endif()

find_path(T_PATH gosam-contrib/libgolem.so PATH_SUFFIXES lib lib64 )
get_filename_component(GOSAM_PREFIX ${T_PATH} DIRECTORY)

# handle the QUIETLY and REQUIRED arguments and set GOSAM_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GoSam DEFAULT_MSG GOSAM_VERSION GOSAM_PREFIX)

mark_as_advanced(GoSam_FOUND)
