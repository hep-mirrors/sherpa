# - Try to find MCFM
# Defines:
#
#  MCFM_FOUND
#  MCFM_PREFIX
if (MCFM_ROOT_DIR OR MCFM_DIR OR (DEFINED ENV{MCFM_ROOT_DIR}) OR (DEFINED ENV{MCFM_DIR}) )
  set(MCFM_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (MCFM_ROOT_DIR)
    list (APPEND MCFM_SEARCH_DIRS "${MCFM_ROOT_DIR}" )
  endif()
  if (MCFM_DIR)
    list (APPEND MCFM_SEARCH_DIRS "${MCFM_DIR}" )
  endif()
  if (DEFINED EVN{MCFM_ROOT_DIR})
    list (APPEND MCFM_SEARCH_DIRS "$ENV{MCFM_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{MCFM_DIR})
    list (APPEND MCFM_SEARCH_DIRS "ENV{MCFM_DIR}" )
  endif()
endif()
if (MCFM_SEARCH_DIRS)
  find_path(MCFM_INCLUDE_DIR  MCFM/CXX_Interface.h PATHS ${MCFM_SEARCH_DIRS} PATH_SUFFIXES include  NO_DEFAULT_PATH)
  find_library(MCFM_LIBRARY NAMES mcfm MCFM PATHS ${MCFM_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_path(MCFM_INCLUDE_DIR  MCFM/CXX_Interface.h PATHS ${MCFM_SEARCH_DIRS} PATH_SUFFIXES include  )
  find_library(MCFM_LIBRARY NAMES mcfm MCFM PATHS_SUFFIXES lib lib64)
endif()
set(MCFM_VERSION Unknown)

get_filename_component(MCFM_PATH ${MCFM_INCLUDE_DIR} DIRECTORY)
# handle the QUIETLY and REQUIRED arguments and set MCFM_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MCFM DEFAULT_MSG MCFM_VERSION  MCFM_LIBRARY MCFM_INCLUDE_DIR MCFM_PATH)

mark_as_advanced(MCFM_FOUND)
