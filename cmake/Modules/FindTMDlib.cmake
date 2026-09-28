# FindTMDlib.cmake
#
# Find the TMDlib library.
#
# This module defines:
#
# TMDlib_FOUND
# TMDlib_INCLUDE_DIRS
# TMDlib_LIBRARIES
# TMDlib_VERSION
#
# Imported target:
#
# TMDlib::tmdlib
#
# User hint:
#
# TMDlib_ROOT_DIR
#
# Expected installation layout:
#
#   TMDlib_ROOT_DIR/
#   ├── bin/
#   │   └── TMDlib-config
#   ├── include/
#   │   └── tmdlib/
#   └── lib/
#       ├── libTMDlib.so
#       └── libapfelxx.so


include(FindPackageHandleStandardArgs)


# ---------------------------------------------------------------------------
# TMDlib root directory
# ---------------------------------------------------------------------------

if(EXISTS "$ENV{TMDlib_ROOT_DIR}")
  file(TO_CMAKE_PATH "$ENV{TMDlib_ROOT_DIR}" TMDlib_ROOT_DIR)
  set(TMDlib_ROOT_DIR
      "${TMDlib_ROOT_DIR}"
      CACHE PATH "Prefix for TMDlib installation."
  )
endif()


# ---------------------------------------------------------------------------
# Try TMDlib-config if no root directory was provided
# ---------------------------------------------------------------------------

if(NOT TMDlib_ROOT_DIR)

  find_program(
    TMDlib_CONFIG_EXECUTABLE
    NAMES TMDlib-config
  )

  if(TMDlib_CONFIG_EXECUTABLE)

    execute_process(
      COMMAND "${TMDlib_CONFIG_EXECUTABLE}" --prefix
      OUTPUT_VARIABLE TMDlib_ROOT_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    set(TMDlib_USE_CONFIG ON)

  endif()

endif()


# ---------------------------------------------------------------------------
# Include directory
# ---------------------------------------------------------------------------

find_path(
  TMDlib_INCLUDE_DIR
  NAMES tmdlib/TMDlib.h
  HINTS
    "${TMDlib_ROOT_DIR}/include"
)


# ---------------------------------------------------------------------------
# TMDlib library
# ---------------------------------------------------------------------------

find_library(
  TMDlib_LIBRARY
  NAMES TMDlib
  HINTS
    "${TMDlib_ROOT_DIR}/lib"
  PATH_SUFFIXES
    Release
    Debug
)


# ---------------------------------------------------------------------------
# APFEL++ library
#
# TMDlib's static library contains unresolved APFEL++ symbols, so APFEL++
# must be available when linking TMDlib.
# ---------------------------------------------------------------------------

find_library(
  TMDlib_APFELXX_LIBRARY
  NAMES apfelxx
  HINTS
    "${TMDlib_ROOT_DIR}/lib"
  PATH_SUFFIXES
    Release
    Debug
)


# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------

if(NOT TMDlib_VERSION)

  find_program(
    TMDlib_CONFIG_EXECUTABLE
    NAMES TMDlib-config
    HINTS
      "${TMDlib_ROOT_DIR}/bin"
  )

  if(TMDlib_CONFIG_EXECUTABLE)

    execute_process(
      COMMAND "${TMDlib_CONFIG_EXECUTABLE}" --version
      OUTPUT_VARIABLE TMDlib_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  endif()

endif()


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

set(TMDlib_INCLUDE_DIRS
    "${TMDlib_INCLUDE_DIR}"
)

set(TMDlib_LIBRARIES
    "${TMDlib_LIBRARY}"
    "${TMDlib_APFELXX_LIBRARY}"
)


# ---------------------------------------------------------------------------
# Standard CMake package handling
# ---------------------------------------------------------------------------

find_package_handle_standard_args(
  TMDlib
  FOUND_VAR TMDlib_FOUND
  REQUIRED_VARS
    TMDlib_INCLUDE_DIR
    TMDlib_LIBRARY
    TMDlib_APFELXX_LIBRARY
  VERSION_VAR
    TMDlib_VERSION
)


# ---------------------------------------------------------------------------
# Imported target
# ---------------------------------------------------------------------------

if(TMDlib_FOUND AND NOT TARGET TMDlib::tmdlib)

  add_library(
    TMDlib::tmdlib
    UNKNOWN
    IMPORTED
  )

  set_target_properties(
    TMDlib::tmdlib
    PROPERTIES
      IMPORTED_LOCATION
        "${TMDlib_LIBRARY}"

      INTERFACE_INCLUDE_DIRECTORIES
        "${TMDlib_INCLUDE_DIRS}"
  )

# TMDlib has unresolved symbols from GSL, LHAPDF, and APFEL++.
#
# GSL and LHAPDF are provided by Sherpa's existing dependency setup.
# APFEL++ is found above in the TMDlib installation itself.

target_link_libraries(TMDlib::tmdlib INTERFACE
  GSL::gsl
  LHAPDF::LHAPDF
  "${TMDlib_APFELXX_LIBRARY}"
)

endif()


# ---------------------------------------------------------------------------
# Advanced cache variables
# ---------------------------------------------------------------------------

mark_as_advanced(
  TMDlib_ROOT_DIR
  TMDlib_CONFIG_EXECUTABLE
  TMDlib_VERSION
  TMDlib_LIBRARY
  TMDlib_APFELXX_LIBRARY
  TMDlib_INCLUDE_DIR
)
