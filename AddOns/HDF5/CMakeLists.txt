# turn off everything we can before fetching HighFive
set(HIGHFIVE_UNIT_TESTS Off CACHE STRING "Enable unit tests (requires Catch2 to be present)")
mark_as_advanced(HIGHFIVE_UNIT_TESTS)
set(HIGHFIVE_EXAMPLES Off)
set(HIGHFIVE_BUILD_DOCS Off)
set(HIGHFIVE_USE_BOOST Off)

# make HighFive available
include(FetchContent)
FetchContent_Declare(
  highfive
  GIT_REPOSITORY https://github.com/highfive-devs/highfive
  GIT_TAG v2.10.1
  # apply patch to fix a CMake deprecation warning for the old minimum
  # version requirement used by HighFive v2.10.1; we simply replace it
  # with our own more recent version requirement
  PATCH_COMMAND git apply "${CMAKE_CURRENT_SOURCE_DIR}/CMake/HDF5-cmake-deprecation.patch"
  # use UPDATE_DISCONNECTED to prevent the patch from being applied (and
  # then fail) again under certain conditions, see
  # https://stackoverflow.com/questions/73718516/modifying-or-injecting-headers-after-fetching-a-git-repo-using-fetchcontent
  UPDATE_DISCONNECTED On
  )
# WORKAROUND: Use the following instead of FetchContent_MakeAvailable,
# because EXCLUDE_FROM_ALL needs to be added via an explicit call to
# add_subdirectory before CMake 3.28 (which we do not require yet).
# See https://stackoverflow.com/questions/65527126/disable-install-for-fetchcontent
FetchContent_GetProperties(highfive)
if(NOT highfive_POPULATED)
  if(POLICY CMP0169)
    # suppress deprecation warning for FetchContent_Populate;
    # CMake 3.28 (see discussion above) is quite new (2023-12-06),
    # so we should wait a bit longer until it's more widely available
    cmake_policy(PUSH)
    cmake_policy(SET CMP0169 OLD)
  endif()
  FetchContent_Populate(highfive)
  if(POLICY CMP0169)
    cmake_policy(POP)
  endif()
  add_subdirectory(${highfive_SOURCE_DIR} ${highfive_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
if (NOT BUILD_SHARED_LIBS)
  install(TARGETS HighFive libheaders libdeps EXPORT PepperTargets)
endif()

set(HDF5_OUTPUT_SOURCES Output_HDF5.C)
set(HDF5_READER_SOURCES Output_HDF5.C HDF5_Reader.C)
sherpa_create_git_info(HDF5 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaHDF5Output SHARED ${HDF5_OUTPUT_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaHDF5Output SherpaHDF5Output)
add_dependencies(SherpaHDF5Output HDF5_Git_Info)
target_link_libraries(SherpaHDF5Output PUBLIC HighFive)
 
add_library(SherpaHDF5Reader SHARED ${HDF5_READER_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaHDF5Reader SherpaHDF5Reader)
add_dependencies(SherpaHDF5Reader HDF5_Git_Info)
target_link_libraries(SherpaHDF5Reader PUBLIC HighFive)
