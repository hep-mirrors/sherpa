# - Try to find RECOLA
# Defines:
#
#  RECOLA_FOUND
#  RECOLA_PREFIX
if (RECOLA_ROOT_DIR OR RECOLA_DIR OR (DEFINED ENV{RECOLA_ROOT_DIR}) OR (DEFINED ENV{RECOLA_DIR}) )
  set(RECOLA_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (RECOLA_ROOT_DIR)
    list (APPEND RECOLA_SEARCH_DIRS "${RECOLA_ROOT_DIR}" )
  endif()
  if (RECOLA_DIR)
    list (APPEND RECOLA_SEARCH_DIRS "${RECOLA_DIR}" )
  endif()
  if (DEFINED EVN{RECOLA_ROOT_DIR})
    list (APPEND RECOLA_SEARCH_DIRS "$ENV{RECOLA_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{RECOLA_DIR})
    list (APPEND RECOLA_SEARCH_DIRS "ENV{RECOLA_DIR}" )
  endif()
endif()
if (RECOLA_SEARCH_DIRS)
  find_path(RECOLA_INCLUDE_DIR recola.h PATHS ${RECOLA_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(RECOLA_LIBRARY NAMES recola PATHS ${RECOLA_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_path(RECOLA_INCLUDE_DIR recola.h PATHS ${RECOLA_SEARCH_DIRS} PATH_SUFFIXES include)
  find_library(RECOLA_LIBRARY NAMES recola PATHS_SUFFIXES lib lib64)
endif()
set(RECOLA_VERSION Unknown)
if (RECOLA_LIBRARY)
  get_filename_component(T_PATH ${RECOLA_LIBRARY} DIRECTORY)
  get_filename_component(RECOLA_PREFIX ${T_PATH} DIRECTORY)
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Recola DEFAULT_MSG RECOLA_VERSION RECOLA_PREFIX RECOLA_LIBRARY RECOLA_INCLUDE_DIR)
mark_as_advanced(Recola_FOUND)
