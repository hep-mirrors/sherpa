
set(SherpaRootTuple_SOURCES Output_RootNtuple.C RootNtuple_Reader.C)
set(SherpaRootTuple_HEADERS Output_RootNtuple.H RootNtuple_Reader.H)
set(SherpaRootOutput_SOURCES Output_RootNtuple.C)
set(SherpaRootOutput_HEADERS Output_RootNtuple.H)
set(SherpaRootInput_SOURCES RootNtuple_Reader.C)
set(SherpaRootInput_HEADERS RootNtuple_Reader.H)
set(SherpaRootEOutput_SOURCES Output_RootNtuple.C)
set(SherpaRootEOutput_HEADERS Output_RootNtuple.H)
set(SherpaRootEInput_SOURCES RootNtuple_Reader.C)
set(SherpaRootEInput_HEADERS RootNtuple_Reader.H)
set(SherpaRootEDOutput_SOURCES Output_RootNtuple.C)
set(SherpaRootEDOutput_HEADERS Output_RootNtuple.H)
set(SherpaRootEDInput_SOURCES RootNtuple_Reader.C)
set(SherpaRootEDInput_HEADERS RootNtuple_Reader.H)

sherpa_create_git_info(SherpaRootTuple ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
list(TRANSFORM SherpaRootOutput_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Root/")
list(APPEND SherpaRootOutput_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

list(TRANSFORM SherpaRootEOutput_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Root/")
list(APPEND SherpaRootEOutput_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

list(TRANSFORM SherpaRootEDOutput_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Root/")
list(APPEND SherpaRootEDOutput_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

list(TRANSFORM SherpaRootInput_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Root/")
list(APPEND SherpaRootInput_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

list(TRANSFORM SherpaRootEInput_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Root/")
list(APPEND SherpaRootEInput_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

list(TRANSFORM SherpaRootEDInput_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Root/")
list(APPEND SherpaRootEDInput_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)

add_library(SherpaRootOutput SHARED ${SherpaRootOutput_SOURCES})
target_link_libraries(SherpaRootOutput PRIVATE ROOT::Tree ROOT::RIO ROOT::Core)
sherpa_handle_shared_library(SherpaRootOutput SherpaRootOutput)

add_library(SherpaERootOutput SHARED ${SherpaRootEOutput_SOURCES})
target_link_libraries(SherpaERootOutput PRIVATE ROOT::Tree ROOT::RIO ROOT::Core)
sherpa_handle_shared_library(SherpaERootOutput SherpaERootOutput)

add_library(SherpaEDRootOutput SHARED ${SherpaRootEDOutput_SOURCES})
target_link_libraries(SherpaEDRootOutput PRIVATE ROOT::Tree ROOT::RIO ROOT::Core)
sherpa_handle_shared_library(SherpaEDRootOutput SherpaEDRootOutput)

add_library(SherpaRootInput SHARED ${SherpaRootInput_SOURCES})
target_link_libraries(SherpaRootInput PRIVATE ROOT::Tree ROOT::RIO ROOT::Core)
sherpa_handle_shared_library(SherpaRootInput SherpaRootInput)

add_library(SherpaERootInput SHARED ${SherpaRootEInput_SOURCES})
target_link_libraries(SherpaERootInput PRIVATE ROOT::Tree ROOT::RIO ROOT::Core)
sherpa_handle_shared_library(SherpaERootInput SherpaERootInput)

add_library(SherpaEDRootInput SHARED ${SherpaRootEDInput_SOURCES})
target_link_libraries(SherpaEDRootInput PRIVATE ROOT::Tree ROOT::RIO ROOT::Core)
sherpa_handle_shared_library(SherpaEDRootInput SherpaEDRootInput)

include(CheckCXXCompilerFlag)
string(REPLACE " " ";" ROOT_CXX_FLAGS_LIST ${ROOT_CXX_FLAGS})
list(REMOVE_ITEM ROOT_CXX_FLAGS_LIST "-stdlib=libc++")
foreach(fl ${ROOT_CXX_FLAGS_LIST})
  CHECK_CXX_COMPILER_FLAG(${fl} COMPILER_SUPPORTS_${fl})
  if(COMPILER_SUPPORTS_${fl})
    target_compile_options(SherpaRootOutput PRIVATE ${fl})
    target_compile_options(SherpaRootInput  PRIVATE ${fl})

    target_compile_options(SherpaERootOutput PRIVATE ${fl})
    target_compile_options(SherpaERootInput  PRIVATE ${fl})

    target_compile_options(SherpaEDRootOutput PRIVATE ${fl})
    target_compile_options(SherpaEDRootInput  PRIVATE ${fl})
  endif()
endforeach()
