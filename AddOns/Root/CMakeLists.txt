#Manually edited
find_package(ROOT REQUIRED) 
set(CMAKE_CXX_FLAGS_BEFORE_ROOT "${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS_BEFORE_ROOT "${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_BEFORE_ROOT "${CMAKE_Fortran_FLAGS}")
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BEFORE_ROOT}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_BEFORE_ROOT}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_BEFORE_ROOT}")
message(STATUS "SHERPA: ROOT_VERSION=${ROOT_VERSION} ROOT_LIBRARIES=${ROOT_LIBRARIES} ROOT_INCLUDE_DIRS=${ROOT_INCLUDE_DIRS} ROOT_CXX_FLAGS=${ROOT_CXX_FLAGS}")

 
set(foo_SOURCES Output_RootNtuple.C RootNtuple_Reader.C)
set(foo_HEADERS Output_RootNtuple.H RootNtuple_Reader.H) 

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


creategitinfo(foo ../Root/)
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
target_link_libraries(SherpaRootOutput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaRootOutput PRIVATE ${ROOT_INCLUDE_DIRS})
sherpa_mpi_link_libraries(SherpaRootOutput)
install(TARGETS SherpaRootOutput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaRootOutput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaRootOutput SOVERSION 0.0.0)
set_target_properties(SherpaRootOutput       PROPERTIES DEFINE_SYMBOL "")

add_library(SherpaERootOutput SHARED ${SherpaRootEOutput_SOURCES})
target_link_libraries(SherpaERootOutput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaERootOutput PRIVATE ${ROOT_INCLUDE_DIRS})
sherpa_mpi_link_libraries(SherpaERootOutput)
install(TARGETS SherpaERootOutput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaERootOutput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaERootOutput SOVERSION 0.0.0)
set_target_properties(SherpaERootOutput       PROPERTIES DEFINE_SYMBOL "")

add_library(SherpaEDRootOutput SHARED ${SherpaRootEDOutput_SOURCES})
target_link_libraries(SherpaEDRootOutput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaEDRootOutput PRIVATE ${ROOT_INCLUDE_DIRS})
sherpa_mpi_link_libraries(SherpaEDRootOutput)
install(TARGETS SherpaEDRootOutput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaEDRootOutput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaEDRootOutput SOVERSION 0.0.0)
set_target_properties(SherpaEDRootOutput       PROPERTIES DEFINE_SYMBOL "")

add_library(SherpaRootInput SHARED ${SherpaRootInput_SOURCES})
target_link_libraries(SherpaRootInput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaRootInput PRIVATE ${ROOT_INCLUDE_DIRS})
sherpa_mpi_link_libraries(SherpaRootInput)
install(TARGETS SherpaRootInput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaRootInput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaRootInput SOVERSION 0.0.0)
set_target_properties(SherpaRootInput       PROPERTIES DEFINE_SYMBOL "")

add_library(SherpaERootInput SHARED ${SherpaRootEInput_SOURCES})
target_link_libraries(SherpaERootInput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaERootInput PRIVATE ${ROOT_INCLUDE_DIRS})
sherpa_mpi_link_libraries(SherpaERootInput)
install(TARGETS SherpaERootInput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaERootInput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaERootInput SOVERSION 0.0.0)
set_target_properties(SherpaERootInput       PROPERTIES DEFINE_SYMBOL "")

add_library(SherpaEDRootInput SHARED ${SherpaRootEDInput_SOURCES})
target_link_libraries(SherpaEDRootInput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaEDRootInput PRIVATE ${ROOT_INCLUDE_DIRS})
sherpa_mpi_link_libraries(SherpaEDRootInput)
install(TARGETS SherpaEDRootInput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaEDRootInput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaEDRootInput SOVERSION 0.0.0)
set_target_properties(SherpaEDRootInput       PROPERTIES DEFINE_SYMBOL "")

#Warning: The ROOT_CXX_FLAGS might be incompatible with the current compiller
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

