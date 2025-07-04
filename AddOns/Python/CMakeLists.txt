find_program(SWIG NAMES swig  PATHS ${SWIG_DIR} PATH_SUFFIXES bin REQUIRED)
if (SWIG)
 message(STATUS "SHERPA: SWIG=${SWIG}")
else()
 message(FATAL_ERROR "SHERPA: SWIG NOT FOUND!")
endif()
set( PYEXTFILES_GEN
    Sherpa.i 
    Exception.i 
    Flavour.i 
    Vec4.i 
    Particle.i 
    Blob.i 
    Blob_List.i 
    MEProcess.i 
    Model_Base.i 
    Random.i 
    Rambo.i 
    Terminator_Objects.i
  )

set(PYTHONSRC)
foreach( file_i ${PYEXTFILES_GEN})
  if ( NOT ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}") )
    add_custom_command(
    PRE_BUILD
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${file_i}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
    )
  endif()
  list( APPEND PYTHONSRC ${CMAKE_CURRENT_BINARY_DIR}/${file_i})
endforeach( file_i )
add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
       DEPENDS ${PYTHONSRC}
       COMMAND ${SWIG} -python -c++  -I${CMAKE_CURRENT_SOURCE_DIR} -I${PROJECT_BINARY_DIR} -I${PROJECT_SOURCE_DIR}   Sherpa.i 
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
Python_add_library(_Sherpa MODULE  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx    )
target_compile_features(_Sherpa PRIVATE cxx_std_11)
set_target_properties(_Sherpa   PROPERTIES      ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               )
if(NOT (CMAKE_SYSTEM_NAME MATCHES Darwin))
  target_link_options(_Sherpa PRIVATE LINKER:--no-as-needed)
endif()
target_include_directories(_Sherpa PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/ ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
install(TARGETS _Sherpa DESTINATION ${SHERPA_Python_SITEARCH}/  COMPONENT python)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Sherpa.py DESTINATION ${SHERPA_Python_SITEARCH}/  COMPONENT python)

set(MEProcess_SOURCES MEProcess.C)
set(MEProcess_HEADERS MEProcess.H)

sherpa_create_git_info(MEProcess ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(MEProcess SHARED ${MEProcess_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(MEProcess MEProcess)

set(_Sherpa_DEPENDENCIES SHERPA::All MEProcess)
if (SHERPA_ENABLE_MPI)
  set(_Sherpa_DEPENDENCIES ${_Sherpa_DEPENDENCIES} MPI::MPI_CXX)
endif ()
target_link_libraries(_Sherpa PRIVATE ${_Sherpa_DEPENDENCIES})
