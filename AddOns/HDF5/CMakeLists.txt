
# make HighFive available
find_package(HighFive QUIET)
if(NOT HighFive_FOUND)
  message(STATUS "HighFive not found. Attempt to install...")

  include(ExternalProject)
  if(HighFive_ARGUMENTS)
    cmake_parse_arguments(HighFive "" "GIT_REPOSITORY;GIT_TAG" ""
      ${HighFive_ARGUMENTS})
  endif()
  if(NOT HighFive_GIT_REPOSITORY)
    set(HighFive_GIT_REPOSITORY https://github.com/highfive-devs/highfive)
  endif()
  if(NOT HighFive_GIT_TAG)
    set(HighFive_GIT_TAG v2.10.1)
  endif()

  set(HighFive_DIR "${CMAKE_INSTALL_PREFIX}/external")
  ExternalProject_Add(
    Lookup-HighFive
    GIT_REPOSITORY ${HighFive_GIT_REPOSITORY}
    GIT_TAG ${HighFive_GIT_TAG}
    PREFIX "${CMAKE_BINARY_DIR}/external"
    #INSTALL_DIR ${EXTERNAL_ROOT}
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${HighFive_DIR}
    -DCMAKE_INSTALL_LIBDIR=${CMAKE_SHARED_LIBRARY_PREFIX}
    -DHIGHFIVE_EXAMPLES=OFF
    -DHIGHFIVE_USE_BOOST=OFF
    -DHIGHFIVE_UNIT_TESTS=OFF
    -DHIGHFIVE_BUILD_DOCS=OFF
    # Wrap download, configure and build steps in a script to log output
    UPDATE_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
  )
  set(HighFive_INCLUDE_DIR "${HighFive_DIR}/include")
  #set(HighFive_LIBRARY_DIR "${HighFive_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}")
  #set(HighFive_LIBRARIES "HighFive")
endif()

sherpa_create_git_info(HDF5 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
set(HDF5_READER_SOURCES HDF5_Reader.C)
add_library(SherpaHDF5Reader SHARED ${HDF5_READER_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaHDF5Reader SherpaHDF5Reader)
if(NOT HighFive_FOUND)
  add_dependencies(SherpaHDF5Reader Lookup-HighFive)
endif()
add_dependencies(SherpaHDF5Reader HDF5_Git_Info)
target_link_libraries(SherpaHDF5Reader PRIVATE ${HDF5_LIBRARIES})
target_include_directories(SherpaHDF5Reader SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS} ${HighFive_INCLUDE_DIR})
