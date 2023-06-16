# - Locate HepMC library
# in a directory defined via  HepMC2_DIR or HEPMC2_DIR.
# Defines:
#
#  HepMC2_FOUND
#  HEPMC2_INCLUDE_DIR
#  HEPMC2_INCLUDE_DIRS (not cached)
#  HEPMC2_LIBRARIES
#  HEPMC2_LIBRARY_DIR (not cached)
if (HepMC2_DIR OR HEPMC2_DIR OR (DEFINED ENV{HepMC2_DIR}) OR (DEFINED ENV{HEPMC2_DIR}) )
  set(HEPMC2_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (HepMC2_DIR)
    list (APPEND HEPMC2_SEARCH_DIRS "${HepMC2_DIR}" )
  endif()
  if (HEPMC2_DIR)
    list (APPEND HEPMC2_SEARCH_DIRS "${HEPMC2_DIR}" )
  endif()
  if (DEFINED ENV{HepMC2_DIR})
    list (APPEND HEPMC2_SEARCH_DIRS "$ENV{HepMC2_DIR}" )
  endif()
  if (DEFINED ENV{HEPMC2_DIR})
    list (APPEND HEPMC2_SEARCH_DIRS "$ENV{HEPMC2_DIR}" )
  endif()
endif()

if (HEPMC2_SEARCH_DIRS)
  find_path(HEPMC2_INCLUDE_DIR HepMC/GenEvent.h PATHS ${HEPMC2_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(HEPMC2_HepMC_LIBRARY NAMES HepMC PATHS ${HEPMC2_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  find_library(HEPMC2_HepMCfio_LIBRARY NAMES HepMCfio  PATHS ${HEPMC2_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_path(HEPMC2_INCLUDE_DIR HepMC/GenEvent.h PATH_SUFFIXES include ../include)
  find_library(HEPMC2_HepMC_LIBRARY NAMES HepMC PATH_SUFFIXES lib lib64 ../lib ../lib64)
  find_library(HEPMC2_HepMCfio_LIBRARY NAMES HepMCfio PATH_SUFFIXES lib lib64 ../lib ../lib64)
endif()

set(HEPMC2_VERSION 0.0.0)
if (HEPMC2_INCLUDE_DIR)
  if (EXISTS ${HEPMC2_INCLUDE_DIR}/HepMC/HepMCDefs.h)
    file(STRINGS ${HEPMC2_INCLUDE_DIR}/HepMC/HepMCDefs.h HEPMC2_VERSION_STRING_CONTENT REGEX "^#define[ ]+HEPMC_VERSION[ ]+\"")
    if (HEPMC2_VERSION_STRING_CONTENT)
      string(REGEX MATCH "[1234567890\.]+[a-zA-Z]*" HEPMC2_VERSION ${HEPMC2_VERSION_STRING_CONTENT})
    endif()
  endif()
endif()

if (HEPMC2_HepMC_LIBRARY)
  set(HEPMC2_LIBRARIES ${HEPMC2_HepMC_LIBRARY})
endif()
if (HEPMC2_HepMCfio_LIBRARY)
  set(HEPMC2_FIO_LIBRARIES ${HEPMC2_HepMCfio_LIBRARY})
  set(HepMC2_FIO_FOUND TRUE)
endif()

set(HEPMC2_INCLUDE_DIRS ${HEPMC2_INCLUDE_DIR})
get_filename_component(HEPMC2_LIBRARY_DIR ${HEPMC2_HepMC_LIBRARY} PATH)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(HepMC2 FOUND_VAR HepMC2_FOUND REQUIRED_VARS HEPMC2_INCLUDE_DIR HEPMC2_LIBRARIES 
                                  VERSION_VAR HEPMC2_VERSION 
                                  HANDLE_VERSION_RANGE
                                  HANDLE_COMPONENTS 
                                 )

mark_as_advanced(HepMC2_FOUND HEPMC2_INCLUDE_DIRS HEPMC2_LIBRARIES HEPMC2_FIO_LIBRARIES)
