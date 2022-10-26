  if (NOT  SHERPA-MC_PYTHON_VERSION)
    find_package (Python  COMPONENTS  Development Interpreter)
  else()
    find_package (Python ${SHERPA-MC_PYTHON_VERSION} EXACT COMPONENTS  Development Interpreter)
  endif()
  if (NOT SHERPA-MC_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
    #set (SHERPA-MC_Python_SITEARCH ${Python_SITEARCH})
    set (SHERPA-MC_Python_SITEARCH ${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/)
    message(STATUS "SHERPA-MC python: SHERPA-MC_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR} defaults to ${SHERPA-MC_Python_SITEARCH}")
  else()
    #set (SHERPA-MC_Python_SITEARCH ${SHERPA-MC_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR}})
    message(STATUS "SHERPA-MC python: SHERPA-MC_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR} is set to to ${SHERPA-MC_Python_SITEARCH}")
  endif()
  if (SHERPA-MC_USE_CYTHON) 
    find_program(CYTHON NAMES cython cython3 cython2)
    message(STATUS "SHERPA-MC python: CYTHON=${CYTHON}")
    if (NOT CYTHON) 
      set(SHERPA-MC_USE_CYTHON OFF)
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
  #list(TRANSFORM PYEXTFILES_GEN PREPEND "${PROJECT_BINARY_DIR}/")

  set( PYEXTFILES

)


  if (NOT SHERPA-MC_USE_CYTHON)
    list (APPEND PYEXTFILES pyext/SHERPA-MC/util.cpp pyext/SHERPA-MC/core.h pyext/SHERPA-MC/core.cpp )

  endif()

  set(PYTHONSRC)

 
  if (SHERPA-MC_USE_CYTHON)

     add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
       DEPENDS ${PYTHONSRC}  ${PYEXTFILES_GEN}
       COMMAND ${CYTHON} ${PYEXTFILES_GEN} --cplus -2  -I ${CMAKE_CURRENT_BINARY_DIR}/pyext -I ${CMAKE_CURRENT_BINARY_DIR}/pyext/SHERPA-MC -I ${CMAKE_CURRENT_BINARY_DIR}/pyext/SHERPA-MC/include  -I ${CMAKE_CURRENT_SOURCE_DIR} -o ${CMAKE_CURRENT_BINARY_DIR}/pyext/SHERPA-MC/core.cpp  
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
  else()
     add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
       DEPENDS ${PYTHONSRC}  ${PYEXTFILES_GEN}
       COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/Sherpa_wrap.cxx ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx
              WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
  
  endif()
  Python_add_library(_Sherpa MODULE  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa_wrap.cxx    )
  set_target_properties(_Sherpa   PROPERTIES      ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/$<0:>
                                               )

  target_link_libraries(_Sherpa PUBLIC ${Python_LIBRARIES} SHERPA )
  target_include_directories(_Sherpa PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/ )
  #add_dependencies(_Sherpa PYTHONSRCTGT)
  install(TARGETS _Sherpa DESTINATION ${SHERPA-MC_Python_SITEARCH}/  COMPONENT python)
  install(FILES Sherpa.py DESTINATION ${SHERPA-MC_Python_SITEARCH}/  COMPONENT python)
  install(FILES ${PROJECT_SOURCE_DIR}/AddOns/DiHiggsNLO/creategrid.py DESTINATION ${SHERPA-MC_Python_SITEARCH}/  COMPONENT python)
