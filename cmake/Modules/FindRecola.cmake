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
  if (DEFINED ENV{RECOLA_ROOT_DIR})
    list (APPEND RECOLA_SEARCH_DIRS "$ENV{RECOLA_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{RECOLA_DIR})
    list (APPEND RECOLA_SEARCH_DIRS "$ENV{RECOLA_DIR}" )
  endif()
endif()
if (RECOLA_SEARCH_DIRS)
  find_path(RECOLA_INCLUDE_DIR recola.h PATHS ${RECOLA_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(RECOLA_LIBRARY NAMES recola PATHS ${RECOLA_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_path(RECOLA_INCLUDE_DIR recola.h PATH_SUFFIXES include ../include)
  find_library(RECOLA_LIBRARY NAMES recola PATH_SUFFIXES lib lib64 ../lib ../lib64)
endif()
if (RECOLA_LIBRARY)
  get_filename_component(T_PATH ${RECOLA_LIBRARY} DIRECTORY)
  get_filename_component(RECOLA_PREFIX ${T_PATH} DIRECTORY)
endif()
find_package(recola CONFIG HINTS ${RECOLA_PREFIX} ${RECOLA_PREFIX}/share QUIET)
if (NOT recola_FOUND)
 set(RECOLA_VERSION 0.0.0)
else()
  set(RECOLA_VERSION ${recola_VERSION})
endif()
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Recola REQUIRED_VARS RECOLA_PREFIX RECOLA_LIBRARY RECOLA_INCLUDE_DIR 
                                 VERSION_VAR RECOLA_VERSION
                                 )

if(Recola_FOUND AND NOT TARGET recola::recola)
    add_library( recola::recola UNKNOWN IMPORTED)
    set_target_properties( recola::recola PROPERTIES
        IMPORTED_LOCATION "${RECOLA_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${RECOLA_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(Recola_FOUND)
