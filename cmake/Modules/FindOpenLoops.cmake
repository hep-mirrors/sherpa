# - Try to find OPENLOOPS
# Defines:
#
#  OPENLOOPS_FOUND
#  OPENLOOPS_INCLUDE_DIR
#  OPENLOOPS_INCLUDE_DIRS (not cached)
#  OPENLOOPS_LIBRARY
#  OPENLOOPS_LIBRARIES (not cached)
#  OPENLOOPS_LIBRARY_DIR (not cached)

if (OPENLOOPS_ROOT_DIR OR OPENLOOPS_DIR OR (DEFINED ENV{OPENLOOPS_ROOT_DIR}) OR (DEFINED ENV{OPENLOOPS_DIR}) )
  set(OPENLOOPS_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (OPENLOOPS_ROOT_DIR)
    list (APPEND OPENLOOPS_SEARCH_DIRS "${OPENLOOPS_ROOT_DIR}" )
  endif()
  if (OPENLOOPS_DIR)
    list (APPEND OPENLOOPS_SEARCH_DIRS "${OPENLOOPS_DIR}" )
  endif()
  if (DEFINED EVN{OPENLOOPS_ROOT_DIR})
    list (APPEND OPENLOOPS_SEARCH_DIRS "$ENV{OPENLOOPS_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{OPENLOOPS_DIR})
    list (APPEND OPENLOOPS_SEARCH_DIRS "ENV{OPENLOOPS_DIR}" )
  endif()
endif()
message(STATUS "${OPENLOOPS_SEARCH_DIRS}")
if (OPENLOOPS_SEARCH_DIRS)
  find_path(OPENLOOPS_PREFIX proclib/channels_public.rinfo PATHS ${OPENLOOPS_SEARCH_DIRS} PATH_SUFFIXES . lib/openloops lib64/openloops  NO_DEFAULT_PATH)
  find_library(OPENLOOPS_LIBRARY NAMES openloops PATHS ${OPENLOOPS_SEARCH_DIRS}  PATH_SUFFIXES . lib lib64 lib/openloops/lib lib64/openloops/lib NO_DEFAULT_PATH)
else()
  find_path(OPENLOOPS_PREFIX proclib/channels_public.rinfo PATH_SUFFIXES . lib/openloops lib64/openloops )
  find_library(OPENLOOPS_LIBRARY NAMES openloops PATHS_SUFFIXES .  lib lib64  lib/openloops/lib lib64/openloops/lib)
endif()

mark_as_advanced(OPENLOOPS_INCLUDE_DIR OPENLOOPS_LIBRARY OPENLOOPS_PREFIX)

# handle the QUIETLY and REQUIRED arguments and set OPENLOOPS_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenLoops DEFAULT_MSG  OPENLOOPS_LIBRARY OPENLOOPS_PREFIX)

set(OPENLOOPS_LIBRARIES ${OPENLOOPS_LIBRARY})
get_filename_component(OPENLOOPS_LIBRARY_DIR ${OPENLOOPS_LIBRARY} PATH)

#set(OPENLOOPS_INCLUDE_DIRS ${OPENLOOPS_INCLUDE_DIR})

mark_as_advanced(OpenLoops_FOUND)
