set(STANDARDENVIRONMENT "SHERPA_CMAKE_PATH=${PROJECT_BINARY_DIR}/outputs/share/SHERPA-MC;LHAPDF_DATA_PATH=${LHAPDF_DATA_PATH};SHERPA_SHARE_PATH=${PROJECT_SOURCE_DIR}/HADRONS++/")
list(APPEND STANDARDENVIRONMENT "LD_LIBRARY_PATH=${LHAPDF_LIBRARY_DIR}:${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC/:$ENV{LD_LIBRARY_PATH}")
list(APPEND STANDARDENVIRONMENT "DYLD_LIBRARY_PATH=${LHAPDF_LIBRARY_DIR}:${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC/:$ENV{DYLD_LIBRARY_PATH}")
if (SHERPA_ENABLE_RIVET)
  list(APPEND STANDARDENVIRONMENT "RIVET_ANALYSIS_PATH=${RIVET_ANALYSIS_PATH};RIVET_DATA_PATH=${RIVET_DATA_PATH}")
  list(APPEND STANDARDENVIRONMENT "PYTHONPATH=${YODA_CONFIG_PYTHONPATH_STRING}:${RIVET_CONFIG_PYTHONPATH_STRING}:$ENV{PYTHONPATH}")
endif()
set(test test-Amegic-LO_Z)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
add_test( NAME Main${test}part1 COMMAND Sherpa ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Sherpa.yaml ME_GENERATORS:Amegic
                              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )
add_test( NAME Main${test}part2 COMMAND ${PROJECT_BINARY_DIR}/AMEGIC++/Main/makelibs
                              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )

add_test( NAME Main${test}part3 COMMAND Sherpa  ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Sherpa.yaml ME_GENERATORS:Amegic
                              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )

SET_TESTS_PROPERTIES( Main${test}part1 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}") 
SET_TESTS_PROPERTIES( Main${test}part2 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")
SET_TESTS_PROPERTIES( Main${test}part3 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}") 

set_tests_properties( Main${test}part2 PROPERTIES DEPENDS Main${test}part1)
set_tests_properties( Main${test}part3 PROPERTIES DEPENDS Main${test}part2)


set(test test-LO_Z)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
add_test( NAME Main${test}part1 COMMAND Sherpa ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Sherpa.yaml 
                              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )
SET_TESTS_PROPERTIES( Main${test}part1 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}") 


if (SHERPA_ENABLE_RIVET)
  set(test test-rivet-LO_Z)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  add_test( NAME Main${test}part1 COMMAND Sherpa ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Sherpa.yaml ANALYSIS=Rivet
                              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )
  SET_TESTS_PROPERTIES( Main${test}part1 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}") 

  add_test( NAME Main${test}part2 COMMAND  ${Python_EXECUTABLE} ${RIVET_MKHTML_EXE} --mc-errs ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Reference.yoda.gz  ${CMAKE_CURRENT_BINARY_DIR}/${test}/Analysis.yoda.gz PLOT:RatioPlotMode=deviation
                                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )
  SET_TESTS_PROPERTIES( Main${test}part2 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}") 

  add_test( NAME Main${test}part3 COMMAND  ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/Examples/CI/plot_deviations ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Reference.yoda.gz ${CMAKE_CURRENT_BINARY_DIR}/${test}/Analysis.yoda.gz
                                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )
  SET_TESTS_PROPERTIES( Main${test}part3 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}") 
  set_tests_properties( Main${test}part2 PROPERTIES DEPENDS Main${test}part1)
  set_tests_properties( Main${test}part3 PROPERTIES DEPENDS Main${test}part2)
endif()
