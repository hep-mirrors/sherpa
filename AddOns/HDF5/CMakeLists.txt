
# make HighFive available
find_package(HighFive QUIET)
if(NOT HighFive_FOUND)
  message(STATUS "HighFive not found. Attempt to fetch...")

  include(FetchContent)

  # make HighFive available
  if(POLICY CMP0135)
    cmake_policy(PUSH)
    cmake_policy(SET CMP0135 NEW)
  endif()
  FetchContent_Declare(
    highfive
    URL https://github.com/highfive-devs/highfive/archive/refs/tags/v3.2.0.tar.gz
    URL_HASH MD5=ebcf6589b3a768696505c5a981f69f58
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/highfive-3.2.0"
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
  if(POLICY CMP0135)
    cmake_policy(POP)
  endif()

endif()

sherpa_create_git_info(HDF5 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
set(HDF5_READER_SOURCES HDF5_Reader.C)
add_library(SherpaHDF5Reader SHARED ${HDF5_READER_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaHDF5Reader SherpaHDF5Reader)
# HighFive is a header-only library. By adding it only
# to the build interface, we can avoid to expose it to consumers of SherpaHDF5Reader,
# and thus do not need to install it. A consequence is that we need
# to explicitly link to HDF5 to make sure that consumers link against it.
target_link_libraries(SherpaHDF5Reader PRIVATE
  $<BUILD_INTERFACE:HighFive::HighFive> HDF5::HDF5)
add_dependencies(SherpaHDF5Reader HDF5_Git_Info)
