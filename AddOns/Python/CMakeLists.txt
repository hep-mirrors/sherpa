if (SHERPA_ENABLE_CYTHON) 
    find_program(CYTHON NAMES cython cython3 cython2 cython-2 cython-3 )
    message(STATUS "SHERPA-MC python: CYTHON=${CYTHON}")
    if (NOT CYTHON) 
      set(SHERPA_ENABLE_CYTHON OFF)
    endif()
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

if (SHERPA_ENABLE_CYTHON)
     add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
       DEPENDS ${PYEXTFILES_GEN}
       COMMAND ${CYTHON} ${PYEXTFILES_GEN} --cplus -2   -I ${CMAKE_CURRENT_SOURCE_DIR} -o ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx 
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
else()
     add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Sherpa_wrap.cxx
       COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/Sherpa_wrap.cxx ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
endif()
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

