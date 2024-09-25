# - Try to find RIVET
# Defines:
#
#  RIVET_FOUND
#  RIVET_INCLUDE_DIR
#  RIVET_INCLUDE_DIRS (not cached)
#  RIVET_LIBRARY
#  RIVET_LIBRARIES (not cached)
#  RIVET_LIBRARY_DIRS (not cached)

if (RIVET_ROOT_DIR OR RIVET_DIR OR (DEFINED ENV{RIVET_ROOT_DIR}) OR (DEFINED ENV{RIVET_DIR}) )
  set(RIVET_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (RIVET_ROOT_DIR)
    list (APPEND RIVET_SEARCH_DIRS "${RIVET_ROOT_DIR}" )
  endif()
  if (RIVET_DIR)
    list (APPEND RIVET_SEARCH_DIRS "${RIVET_DIR}" )
  endif()
  if (DEFINED ENV{RIVET_ROOT_DIR})
    list (APPEND RIVET_SEARCH_DIRS "$ENV{RIVET_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{RIVET_DIR})
    list (APPEND RIVET_SEARCH_DIRS "$ENV{RIVET_DIR}" )
  endif()
endif()

if (RIVET_SEARCH_DIRS)
  find_program(RIVET_MKHTML_EXE NAMES rivet-mkhtml PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES bin NO_DEFAULT_PATH )
  find_program(RIVET_CONFIG_EXE NAMES rivet-config PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES bin NO_DEFAULT_PATH )
  find_program(RIVET_EXE NAMES rivet PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES bin NO_DEFAULT_PATH )
  find_path(RIVET_DATA_PATH NAMES ATLAS_2012_I1094061.yoda.gz ATLAS_2012_I1094061.yoda PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES share/Rivet/ NO_DEFAULT_PATH)
  find_path(RIVET_INCLUDE_DIR Rivet/Rivet.hh PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(RIVET_LIBRARY NAMES Rivet PATHS ${RIVET_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_program(RIVET_MKHTML_EXE NAMES rivet-mkhtml)
  find_program(RIVET_CONFIG_EXE NAMES rivet-config)
  find_program(RIVET_EXE NAMES rivet)
  find_path(RIVET_DATA_PATH NAMES ATLAS_2012_I1094061.yoda.gz ATLAS_2012_I1094061.yoda PATH_SUFFIXES share/Rivet/ ../share/Rivet/)
  find_path(RIVET_INCLUDE_DIR Rivet/Rivet.hh PATH_SUFFIXES include ../include)
  find_library(RIVET_LIBRARY NAMES Rivet PATH_SUFFIXES lib lib64 ../lib ../lib64)
endif()
if ( NOT RIVET_DATA_PATH AND EXISTS ${RIVET_INCLUDE_DIR}/../share/Rivet)
  set(RIVET_DATA_PATH ${RIVET_INCLUDE_DIR}/../share/Rivet)
endif()
set(RIVET_VERSION 0.0.0)
if (RIVET_INCLUDE_DIR)
  if (EXISTS ${RIVET_INCLUDE_DIR}/Rivet/Config/RivetConfig.hh)
    file(STRINGS ${RIVET_INCLUDE_DIR}/Rivet/Config/RivetConfig.hh RIVET_VERSION_STRING_CONTENT REGEX "^#define[ ]+RIVET_VERSION[ ]+\"" )
    if (RIVET_VERSION_STRING_CONTENT)
      string(REGEX MATCH "[1234567890.]+[a-zA-Z]*" RIVET_VERSION ${RIVET_VERSION_STRING_CONTENT})
      # TODO:
      # to be done properly (use Rivet's RIVET_VERSION_CODE instead)
      string(REPLACE "."  ";" RIVET_VERSION_LIST "${RIVET_VERSION}")
      list(GET RIVET_VERSION_LIST 0 RIVET_VERSION_MAJOR)
      list(GET RIVET_VERSION_LIST 1 RIVET_VERSION_MINOR)
      list(GET RIVET_VERSION_LIST 2 RIVET_VERSION_PATCH)
    endif()
    file(STRINGS ${RIVET_INCLUDE_DIR}/Rivet/Config/RivetConfig.hh RIVET_HEPMC_VERSION_STRING REGEX "^#define[ ]+RIVET_ENABLE_HEPMC_3[ ]+" )
    if (RIVET_HEPMC_VERSION_STRING)
      set(Rivet_HEPMC2_FOUND FALSE)
      set(Rivet_HEPMC3_FOUND TRUE)
    else()
      set(Rivet_HEPMC2_FOUND TRUE)
      set(Rivet_HEPMC3_FOUND FALSE)
    endif()
  endif()
endif()

set(RIVET_CONFIG_LIBS_STRING)
set(RIVET_CONFIG_CPPFLAGS_STRING)
set(RIVET_CONFIG_CPPFLAGS_DIRS "")
set(RIVET_CONFIG_LIBS "")
set( RIVET_CONFIG_LIB_DIRS "")
if (RIVET_CONFIG_EXE)
  execute_process(COMMAND ${RIVET_CONFIG_EXE} --libs
                  OUTPUT_VARIABLE RIVET_CONFIG_LIBS_STRING
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REPLACE " " ";" TEMP_RIVET_LIBS  ${RIVET_CONFIG_LIBS_STRING})
  foreach (fl ${TEMP_RIVET_LIBS})
  if ("${fl}" MATCHES "-l.*")
    string(REPLACE "-l" ""  flx "${fl}")
    list(APPEND RIVET_CONFIG_LIBS "${flx}")
  endif()
  if ("${fl}" MATCHES "-L.*")
    string(REPLACE "-L" ""  fly "${fl}")
    list(APPEND RIVET_CONFIG_LIB_DIRS "${fly}")
  endif()
  endforeach()
  execute_process(COMMAND ${RIVET_CONFIG_EXE} --cppflags
                  OUTPUT_VARIABLE RIVET_CONFIG_CPPFLAGS_STRING
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REPLACE "-I" "" TEMP_RIVET_CONFIG_CPPFLAGS_DIRS  ${RIVET_CONFIG_CPPFLAGS_STRING})
  string(REPLACE " " ";" RIVET_CONFIG_CPPFLAGS_DIRS_LIST  ${TEMP_RIVET_CONFIG_CPPFLAGS_DIRS})
  set(RIVET_CONFIG_CPPFLAGS_DIRS )
  foreach( X IN LISTS RIVET_CONFIG_CPPFLAGS_DIRS_LIST)
    if (EXISTS ${X})
       LIST(APPEND RIVET_CONFIG_CPPFLAGS_DIRS ${X})
    endif()
  endforeach()
  execute_process(COMMAND ${RIVET_CONFIG_EXE} --pythonpath
                  OUTPUT_VARIABLE RIVET_CONFIG_PYTHONPATH_STRING
                  OUTPUT_STRIP_TRAILING_WHITESPACE) 
endif()

if (RIVET_VERSION VERSION_GREATER_EQUAL 4.0.0)
    set(RIVET_YODA_MIN_VERSION "2.0.0")
    set(RIVET_HEPMC3_MIN_VERSION "3.2.6")    
    set(RIVET_MKHTML_ARGS )
else()
    set(RIVET_YODA_MIN_VERSION "1.8.0")
    set(RIVET_HEPMC3_MIN_VERSION "3.2.3")
    set(RIVET_MKHTML_ARGS "--mc-errs")
endif()
message(STATUS "SHERPA: WARNING: YODA with RIVET_YODA_MIN_VERSION=${RIVET_YODA_MIN_VERSION} will be requested by Rivet.")
find_package(YODA ${RIVET_YODA_MIN_VERSION} REQUIRED)

mark_as_advanced(RIVET_INCLUDE_DIR RIVET_LIBRARY RIVET_EXE RIVET_CONFIG_LIBS_STRING 
                               RIVET_CONFIG_CPPFLAGS_STRING 
                               RIVET_CONFIG_CPPFLAGS_DIRS
                               RIVET_CONFIG_LIBS
                               RIVET_CONFIG_LIB_DIRS
                               RIVET_MKHTML_ARGS
                               )
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Rivet HANDLE_COMPONENTS REQUIRED_VARS RIVET_INCLUDE_DIR RIVET_LIBRARY 
                               RIVET_DATA_PATH
                               RIVET_CONFIG_LIBS_STRING 
                               RIVET_CONFIG_CPPFLAGS_STRING 
                               RIVET_CONFIG_CPPFLAGS_DIRS
                               RIVET_CONFIG_LIBS
                               RIVET_CONFIG_LIB_DIRS
                               RIVET_YODA_MIN_VERSION
                               RIVET_HEPMC3_MIN_VERSION
                               VERSION_VAR RIVET_VERSION
                               )

set(RIVET_LIBRARIES ${RIVET_LIBRARY})
get_filename_component(RIVET_LIBRARY_DIRS ${RIVET_LIBRARY} PATH)
get_filename_component(RIVET_ANALYSIS_PATH ${RIVET_LIBRARY} PATH)
set(RIVET_ANALYSIS_PATH ${RIVET_ANALYSIS_PATH}/Rivet)
set(RIVET_INCLUDE_DIRS ${RIVET_INCLUDE_DIR})



if(RIVET_FOUND AND NOT TARGET Rivet::Rivet)
    add_library(Rivet::Rivet UNKNOWN IMPORTED)
    if (RIVET_VERSION VERSION_GREATER_EQUAL 4.0.0)
      set_target_properties(Rivet::Rivet PROPERTIES
        IMPORTED_LOCATION "${RIVET_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${RIVET_INCLUDE_DIR};${RIVET_CONFIG_CPPFLAGS_DIRS}"
        INTERFACE_COMPILE_FEATURES cxx_std_17
      )
  else()
      set_target_properties(Rivet::Rivet PROPERTIES
        IMPORTED_LOCATION "${RIVET_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${RIVET_INCLUDE_DIR};${RIVET_CONFIG_CPPFLAGS_DIRS}"
        INTERFACE_COMPILE_FEATURES cxx_std_14
      )
  endif()
endif()

mark_as_advanced(RIVET_FOUND)
