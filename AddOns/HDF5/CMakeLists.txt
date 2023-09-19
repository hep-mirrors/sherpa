cmake_minimum_required(VERSION 3.20)
find_package(HDF5 REQUIRED)

set(SherpaHDF5_SOURCES Output_HDF5.C )
sherpa_create_git_info(SherpaHDF5)

add_library(SherpaHDF5 SHARED ${SherpaHDF5_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaHDF5 SherpaHDF5)
target_include_directories(SherpaHDF5 PRIVATE ${HDF5_INCLUDE_DIR})
target_link_libraries(SherpaHDF5 PRIVATE HDF5::HDF5)
