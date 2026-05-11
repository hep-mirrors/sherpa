
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

# Reader
set(HDF5_READER_SOURCES HDF5_Reader.C)
add_library(SherpaHDF5Reader SHARED ${HDF5_READER_SOURCES}
  ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaHDF5Reader SherpaHDF5Reader)
# HighFive is a header-only library. By adding it only
# to the build interface, we can avoid to expose it to consumers of SherpaHDF5Reader,
# and thus do not need to install it. A consequence is that we need
# to explicitly link to HDF5 to make sure that consumers link against it.
target_link_libraries(SherpaHDF5Reader PRIVATE
  $<BUILD_INTERFACE:HighFive::HighFive> HDF5::HDF5)
add_dependencies(SherpaHDF5Reader HDF5_Git_Info)

# Output
set(HDF5_OUTPUT_SOURCES Output_HDF5.C)
add_library(SherpaHDF5Output SHARED ${HDF5_OUTPUT_SOURCES}
  ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaHDF5Output SherpaHDF5Output)
# HighFive is a header-only library. By adding it only
# to the build interface, we can avoid to expose it to consumers of SherpaHDF5Output,
# and thus do not need to install it. A consequence is that we need
# to explicitly link to HDF5 to make sure that consumers link against it.
target_link_libraries(SherpaHDF5Output PRIVATE
  $<BUILD_INTERFACE:HighFive::HighFive> HDF5::HDF5)
add_dependencies(SherpaHDF5Output HDF5_Git_Info)

# Pepper reader. Reads parton-level events from in-memory HDF5
# databases that are filled by the Pepper library.
if (SHERPA_ENABLE_PEPPER)
  set(PEPPER_READER_SOURCES Pepper_Reader.C)
  add_library(SherpaPepperReader SHARED ${PEPPER_READER_SOURCES}
    ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
  sherpa_handle_shared_library(SherpaPepperReader SherpaPepperReader)
  # HighFive is header-only (build-interface only, same rationale as
  # for SherpaHDF5Reader); Pepper is linked publicly because its API
  # is part of the reader's translation unit.
  target_link_libraries(SherpaPepperReader PRIVATE
    $<BUILD_INTERFACE:HighFive::HighFive> HDF5::HDF5 Pepper::pepperlib)
  add_dependencies(SherpaPepperReader HDF5_Git_Info)
endif()
