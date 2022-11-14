# Find libzip library and headers
#
# The module defines the following variables:
#
# ::
#
#   LibZip_FOUND             - true if libzip was found
#   LibZip_INCLUDE_DIRS      - include search path
#   LibZip_LIBRARIES         - libraries to link
#   LibZip_VERSION           - libzip 3-component version number

find_package(PkgConfig)
pkg_check_modules(PC_LIBZIP QUIET libzip)

if (PC_LIBZIP_VERSION)
 set(LibZip_VERSION ${PC_LIBZIP_VERSION})
else()
 set(LibZip_VERSION Unknown)
endif()

find_path(LibZip_INCLUDE_DIR zip.h
  HINTS ${PC_LIBZIP_INCLUDEDIR})

# Contains the version of libzip:
find_path(LibZip_INCLUDE_CONF_DIR zipconf.h
  HINTS ${PC_LIBZIP_INCLUDE_DIRS})

find_library(LibZip_LIBRARIES
  NAMES zip libzip
  HINTS ${PC_LIBZIP_LIBDIR})

set(LibZip_INCLUDE_DIRS ${LibZip_INCLUDE_DIR} ${LibZip_INCLUDE_CONF_DIR})

if ( NOT LibZip_INCLUDE_DIR OR NOT LibZip_INCLUDE_CONF_DIR OR NOT LibZip_LIBRARIES)
  if (LibZip_ROOT_DIR OR LibZip_DIR OR (DEFINED ENV{LibZip_ROOT_DIR}) OR (DEFINED ENV{LibZip_DIR}) )
    set(LibZip_SEARCH_DIRS "" CACHE STRING "" FORCE)
    if (LibZip_ROOT_DIR)
      list (APPEND LibZip_SEARCH_DIRS "${LibZip_ROOT_DIR}" )
    endif()
    if (LibZip_DIR)
      list (APPEND LibZip_SEARCH_DIRS "${LibZip_DIR}" )
    endif()
    if (DEFINED EVN{LibZip_ROOT_DIR})
      list (APPEND LibZip_SEARCH_DIRS "$ENV{LibZip_ROOT_DIR}" )
    endif()
    if (DEFINED ENV{LibZip_DIR})
      list (APPEND LibZip_SEARCH_DIRS "ENV{LibZip_DIR}" )
    endif()
  endif()
  if (LibZip_SEARCH_DIRS)
    find_path(LibZip_INCLUDE_DIR zip.h PATHS ${LibZip_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
    find_path(LibZip_INCLUDE_CONF_DIR zipconf.h PATHS ${LibZip_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
    find_library(LibZip_LIBRARIES NAMES zip libzip PATHS ${LibZip_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  else()
    find_path(LibZip_INCLUDE_DIR zip.h PATHS ${LibZip_SEARCH_DIRS} PATH_SUFFIXES include)
    find_path(LibZip_INCLUDE_CONF_DIR zipconf.h PATHS ${LibZip_SEARCH_DIRS} PATH_SUFFIXES include)
    find_library(LibZip_LIBRARIES NAMES zip libzip PATHS_SUFFIXES lib lib64)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibZip
                                  FOUND_VAR LibZip_FOUND
                                  REQUIRED_VARS LibZip_LIBRARIES LibZip_INCLUDE_DIR LibZip_INCLUDE_CONF_DIR
                                  VERSION_VAR LibZip_VERSION)

if(LibZip_FOUND AND NOT TARGET LibZip::LibZip)
    add_library(LibZip::LibZip UNKNOWN IMPORTED)
    set_target_properties(LibZip::LibZip PROPERTIES
        IMPORTED_LOCATION "${LibZip_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${LibZip_INCLUDE_DIRS}"
    )
endif()

mark_as_advanced(LibZip_INCLUDE_DIR LibZip_INCLUDE_CONF_DIR)
