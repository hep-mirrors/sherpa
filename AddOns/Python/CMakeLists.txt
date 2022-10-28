
    find_program(SWIG NAMES swig  REQUIRED)
    message(STATUS "SHERPA-MC python: SWIG=${SWIG}")
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

#if (SHERPA_ENABLE_SWIG)
  set(PYTHONSRC)
  foreach( file_i ${PYEXTFILES_GEN})
    add_custom_command(
    PRE_BUILD
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${file_i}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
    )
    list( APPEND PYTHONSRC ${CMAKE_CURRENT_BINARY_DIR}/${file_i})
  endforeach( file_i )
     add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
       DEPENDS ${PYTHONSRC}
       COMMAND ${SWIG} -c++ -python -I${CMAKE_CURRENT_SOURCE_DIR} -I${PROJECT_BINARY_DIR} -I${PROJECT_SOURCE_DIR}   ${PYEXTFILES_GEN}
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
#else()
#     add_custom_command(
#       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
#       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Sherpa_wrap.cxx
#       COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/Sherpa_wrap.cxx ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
#       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
#       VERBATIM)
#endif()
Python_add_library(_Sherpa MODULE  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx    )
set_target_properties(_Sherpa   PROPERTIES      ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               )
target_include_directories(_Sherpa PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/ )
install(TARGETS _Sherpa DESTINATION ${SHERPA-MC_Python_SITEARCH}/  COMPONENT python)
install(FILES Sherpa.py DESTINATION ${SHERPA-MC_Python_SITEARCH}/  COMPONENT python)
install(FILES ${PROJECT_SOURCE_DIR}/AddOns/DiHiggsNLO/creategrid.py DESTINATION ${SHERPA-MC_Python_SITEARCH}/  COMPONENT python)

add_library(MEProcess SHARED MEProcess.C)
install(TARGETS MEProcess DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MEProcess       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MEProcess SOVERSION 0.0.0)

set(PYLIBDIR ${SHERPA-MC_Python_SITEARCH})
configure_file(${PROJECT_SOURCE_DIR}//Examples/API/Events/test.py.in ${PROJECT_SOURCE_DIR}//Examples/API/Events/test.py  @ONLY)
configure_file(${PROJECT_SOURCE_DIR}//Examples/API/ME2-Python/test.py.in ${PROJECT_SOURCE_DIR}//Examples/API/ME2-Python/test.py @ONLY)
configure_file(${PROJECT_SOURCE_DIR}//Examples/API/MPIEvents/test.py.in ${PROJECT_SOURCE_DIR}//Examples/API/MPIEvents/test.py @ONLY)

