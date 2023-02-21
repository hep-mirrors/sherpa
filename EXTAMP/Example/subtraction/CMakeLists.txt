project(ExampleME LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.11.0)
find_package(SHERPA-MC REQUIRED)
add_library(ExampleME SHARED ExampleME.C ExampleME.H CPPProcess_born.h  CPPProcess_real.h  HelAmps_sm.h  Parameters_sm.h  rambo.h  read_slha.h 
                                                     CPPProcess_born.cc  CPPProcess_real.cc  HelAmps_sm.cc  Parameters_sm.cc  rambo.cc  read_slha.cc)
amegic_handle_shared_library(ExampleME ExampleME)
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
install(TARGETS ExampleME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

