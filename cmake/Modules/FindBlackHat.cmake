# - Try to find BLACKHAT
# Defines:
#
#  BLACKHAT_FOUND
#  BLACKHAT_INCLUDE_DIR
#  BLACKHAT_INCLUDE_DIRS (not cached)
#  BLACKHAT_LIBRARY
#  BLACKHAT_LIBRARIES (not cached)
#  BLACKHAT_LIBRARY_DIR (not cached)

if (BLACKHAT_ROOT_DIR OR BLACKHAT_DIR OR (DEFINED ENV{BLACKHAT_ROOT_DIR}) OR (DEFINED ENV{BLACKHAT_DIR}) )
  set(BLACKHAT_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (BLACKHAT_ROOT_DIR)
    list (APPEND BLACKHAT_SEARCH_DIRS "${BLACKHAT_ROOT_DIR}" )
  endif()
  if (BLACKHAT_DIR)
    list (APPEND BLACKHAT_SEARCH_DIRS "${BLACKHAT_DIR}" )
  endif()
  if (DEFINED ENV{BLACKHAT_ROOT_DIR})
    list (APPEND BLACKHAT_SEARCH_DIRS "$ENV{BLACKHAT_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{BLACKHAT_DIR})
    list (APPEND BLACKHAT_SEARCH_DIRS "$ENV{BLACKHAT_DIR}" )
  endif()
endif()
if (BLACKHAT_SEARCH_DIRS)
  find_path(BLACKHAT_INCLUDE_DIR blackhat/BH_interface.h PATHS ${BLACKHAT_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(BLACKHAT_LIBRARY NAMES BH PATHS ${BLACKHAT_SEARCH_DIRS}  PATH_SUFFIXES lib/blackhat lib64/blackhat NO_DEFAULT_PATH)
else()
  find_path(BLACKHAT_INCLUDE_DIR blackhat/BH_interface.h PATH_SUFFIXES include ../include)
  find_library(BLACKHAT_LIBRARY NAMES BH PATH_SUFFIXES lib/blackhat lib64/blackhat ../lib/blackhat ../lib64/blackhat)
endif()

mark_as_advanced(BLACKHAT_INCLUDE_DIR BLACKHAT_LIBRARY)
get_filename_component(BLACKHAT_PATH ${BLACKHAT_INCLUDE_DIR} DIRECTORY)

find_program(BHCONFIG NAMES blackhat-config PATH ${BLACKHAT_PATH}/bin)
execute_process(COMMAND ${BHCONFIG} --version
                OUTPUT_VARIABLE BLACKHAT_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT BLACKHAT_VERSION)
  set(BLACKHAT_VERSION Unknown)
endif()

set(BLACKHAT_LIBRARIES )
get_filename_component(BLACKHAT_LIBRARY_DIR ${BLACKHAT_LIBRARY} PATH)

set (ALLBHTOFIND BH Interface assembly OLA CutPart RatPart ratext Ampl_eval Rateval Cuteval Cut_wCI BG BHcore Integrals Spinors)
set (BLACKHAT_LIBRARIES_ALL TRUE)
foreach(fl IN LISTS ALLBHTOFIND)
  find_library(BLACKHAT_LIBRARY_${fl} NAMES ${fl} PATHS ${BLACKHAT_LIBRARY_DIR}  NO_DEFAULT_PATH)
  if (BLACKHAT_LIBRARY_${fl})
    list(APPEND BLACKHAT_LIBRARIES ${BLACKHAT_LIBRARY_${fl}})
  else()
    set(BLACKHAT_LIBRARIES_ALL FALSE)
  endif()
endforeach()


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BlackHat REQUIRED_VARS BLACKHAT_INCLUDE_DIR BLACKHAT_LIBRARIES BLACKHAT_PATH BLACKHAT_LIBRARIES_ALL
                                  VERSION_VAR BLACKHAT_VERSION 
                                 )
if(BlackHat_FOUND)
add_library(BlackHat::BlackHatAll INTERFACE IMPORTED)
  foreach(fl IN LISTS ALLBHTOFIND)
    add_library( BlackHat::${fl} UNKNOWN IMPORTED)
    set_target_properties( BlackHat::${fl} PROPERTIES
        IMPORTED_LOCATION "${BLACKHAT_LIBRARY_${fl}}"
        INTERFACE_INCLUDE_DIRECTORIES "${BLACKHAT_INCLUDE_DIR}"
    )
  target_link_libraries( BlackHat::BlackHatAll INTERFACE BlackHat::${fl})
  endforeach()
endif()

set(BLACKHAT_INCLUDE_DIRS ${BLACKHAT_INCLUDE_DIR})

mark_as_advanced(BlackHat_FOUND)
