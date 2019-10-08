cmake_minimum_required (VERSION 2.6)

# project name and version
project(test_recola)
set( version 1.0.0 )

if(RECOLA_PREFIX)
  list(APPEND CMAKE_PREFIX_PATH ${RECOLA_PREFIX})
endif()

find_package(Yaml-cpp REQUIRED)
if(${Yaml-cpp_FOUND})
  message(STATUS "Yaml-correctly found")
  set(YAML_PREFIX "${YAML_CPP_INCLUDE_DIR}/..")
endif()
if(NOT YAML_PREFIX)
  message( FATAL_ERROR
    "yaml-cpp path not found automatically. If you have it installed, set it via -DYAML_PREFIX=<path-to-yaml-cpp>.
If you don't have yaml-cpp installed, please install it and try again" )
endif()

find_package(recola REQUIRED)
if(${recola_FOUND})
  message(STATUS "Recola-correctly found")
  set(RECOLA_PREFIX "${recola_INCLUDE_DIR}/..")
endif()
if(NOT RECOLA_PREFIX)
  message( FATAL_ERROR
    "recola path not found automatically. If you have it installed, set it via -DRECOLA_PREFIX=<path-to-recola-cpp>.
If you don't have recola-cpp installed, please install it and try again" )
endif()

set(CMAKE_CXX_FLAGS "-g -O3 -fPIC -Wall")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR} source)
include_directories ("${YAML_PREFIX}/include")
include_directories ("${RECOLA_PREFIX}/include")

link_directories ("${YAML_PREFIX}/lib")
link_directories ("${RECOLA_PREFIX}")
link_libraries(yaml-cpp recola)


add_executable(test-recola ${source})
