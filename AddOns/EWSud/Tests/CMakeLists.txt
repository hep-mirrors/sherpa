set( EWSud_TESTS  eePP eeWW eeZP eeZZ eebb eedd eemumu eett eeuu duWZ duWP uuPG udWG uuZG)
foreach ( test ${EWSud_TESTS} )
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  add_test( NAME EWSud${test} COMMAND Sherpa  "RUNDATA: [ ${CMAKE_CURRENT_SOURCE_DIR}/Setups/Base.yaml, ${CMAKE_CURRENT_SOURCE_DIR}/Setups/${test}.yaml]"
                              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
           )
  SET_TESTS_PROPERTIES( EWSud${test} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")
endforeach ( test ${EWSud_TESTS} )
