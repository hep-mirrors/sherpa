# - Try to find HZTOOL
# Defines:
#
#  HZTOOL_FOUND
#  HZTOOL_PREFIX
if (HZTOOL_ROOT_DIR OR HZTOOL_DIR OR (DEFINED ENV{HZTOOL_ROOT_DIR}) OR (DEFINED ENV{HZTOOL_DIR}) )
  set(HZTOOL_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (HZTOOL_ROOT_DIR)
    list (APPEND HZTOOL_SEARCH_DIRS "${HZTOOL_ROOT_DIR}" )
  endif()
  if (HZTOOL_DIR)
    list (APPEND HZTOOL_SEARCH_DIRS "${HZTOOL_DIR}" )
  endif()
  if (DEFINED EVN{HZTOOL_ROOT_DIR})
    list (APPEND HZTOOL_SEARCH_DIRS "$ENV{HZTOOL_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{HZTOOL_DIR})
    list (APPEND HZTOOL_SEARCH_DIRS "ENV{HZTOOL_DIR}" )
  endif()
endif()
if (HZTOOL_SEARCH_DIRS)
  find_path(HZTOOL_INCLUDE_DIR  heracmn.inc PATHS ${HZTOOL_SEARCH_DIRS} PATH_SUFFIXES include/hztool  NO_DEFAULT_PATH)
  find_library(HZTOOL_LIBRARY NAMES hztool PATHS ${HZTOOL_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
find_path(HZTOOL_INCLUDE_DIR  heracmn.inc PATH_SUFFIXES include/hztool  ../include/hztool)
  find_library(HZTOOL_LIBRARY NAMES hztool PATH_SUFFIXES lib lib64 ../lib ../lib64)
endif()
set(HZTOOL_VERSION 0.0.0)
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HZTool REQUIRED_VARS HZTOOL_INCLUDE_DIR HZTOOL_LIBRARY
                                 VERSION_VAR HZTOOL_VERSION
                                 )

mark_as_advanced(HZTool_FOUND)
