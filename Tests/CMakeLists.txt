if (SHERPA_ENABLE_RIVET)
  list(APPEND STANDARDENVIRONMENT "RIVET_ANALYSIS_PATH=${RIVET_ANALYSIS_PATH};RIVET_DATA_PATH=${RIVET_DATA_PATH}")
  list(APPEND STANDARDENVIRONMENT "PYTHONPATH=${YODA_CONFIG_PYTHONPATH_STRING}:${RIVET_CONFIG_PYTHONPATH_STRING}:$ENV{PYTHONPATH}")
endif()
if (SHERPA_ENABLE_INTEGRATION_TESTS)
  set(test test-Amegic-LO_Z)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  add_test( NAME Main${test}part1 COMMAND Sherpa ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Sherpa.yaml ME_GENERATORS:Amegic
                                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
             )
  add_test( NAME Main${test}part2 COMMAND ${makelibsscript}
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

  if (SHERPA_ENABLE_PYTHIA8)
    set(test test-LO_Z-Pythia8)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
    add_test( NAME Main${test}part1 COMMAND Sherpa  "FRAGMENTATION: Pythia8"  ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Sherpa.yaml 
                                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
               )
    SET_TESTS_PROPERTIES( Main${test}part1 PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT};PYTHIA8DATA=${PYTHIA8_XMLDOC_DIR}")
  endif()

  # ==============================
  # Declare test "test-rivet-LO_Z"
  # ==============================
  #
  # Run a simple Z production set-up and compare its Rivet analysis results
  # against previously generated reference results.
  #
  if(SHERPA_ENABLE_RIVET)
    # create test directory within the build directory
    set(test test-rivet-LO_Z)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})

    # run Sherpa (part 1/3)
    add_test(
      NAME Main${test}part1
      COMMAND Sherpa ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Sherpa.yaml ANALYSIS=Rivet
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
      )
    set_tests_properties(Main${test}part1
      PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")

    # run rivet-mkhtml (part 2/3)
    add_test(
      NAME Main${test}part2
      COMMAND ${Python_EXECUTABLE} ${RIVET_MKHTML_EXE} ${RIVET_MKHTML_ARGS}
        ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Reference.yoda.gz
        ${CMAKE_CURRENT_BINARY_DIR}/${test}/Analysis.yoda.gz
        PLOT:RatioPlotMode=deviation
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
      )
    set_tests_properties(Main${test}part2 PROPERTIES
      ENVIRONMENT "${STANDARDENVIRONMENT}"
      DEPENDS Main${test}part1)

    # plot and check deviations (part 3/3)
    add_test(
      NAME Main${test}part3
      COMMAND ${Python_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/Examples/CI/plot_deviations
        ${PROJECT_SOURCE_DIR}/Examples/CI/LO_Z/Reference.yoda.gz
        ${CMAKE_CURRENT_BINARY_DIR}/${test}/Analysis.yoda.gz
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
      )
    set_tests_properties(Main${test}part3 PROPERTIES
      ENVIRONMENT "${STANDARDENVIRONMENT}"
      DEPENDS Main${test}part2)
  endif()


  # ================================================
  # Declare tests "test-rivet-[LO|NLO]_Z_Variations"
  # ================================================
  #
  # Run a simple Z production set-up with on-the-fly variations and compare its
  # Rivet analysis results with runs that use explicit variations, both at LO and
  # at NLO.
  #
  if(SHERPA_ENABLE_RIVET)
    # macro that declares the tests, order is either LO or NLO
    macro(declare_variations_test order)

      # create test directory within the build directory
      set(test test-rivet-${order}_Z_Variations)
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})

      # do initialization and integration run (part A)
      # doing this separately ensures that all runs below start with the same
      # Results.db and random state
      add_test(
        NAME Main${test}partA
        COMMAND Sherpa -e 0 ${PROJECT_SOURCE_DIR}/Examples/CI/FO_Z_Variations/Base.yaml ${PROJECT_SOURCE_DIR}/Examples/CI/FO_Z_Variations/${order}.yaml
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
        )

      # do on-the-fly run
      add_test(
        NAME Main${test}partB
        COMMAND Sherpa -A OTF ${PROJECT_SOURCE_DIR}/Examples/CI/FO_Z_Variations/Base.yaml ${PROJECT_SOURCE_DIR}/Examples/CI/FO_Z_Variations/${order}.yaml
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
        )
      set_tests_properties(Main${test}partA Main${test}partB
        PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")
      set_tests_properties(Main${test}partB PROPERTIES DEPENDS Main${test}partA)

      # now do explicit variation runs
      macro(run_var NN a b c d)
        add_test(
          NAME Main${test}partC${NN}
          COMMAND Sherpa SCALE_VARIATIONS:None RSF:=${a} FSF:=${b} -A Explicit__MUR=${c}__MUF=${d}__LHAPDF=93300
            ${PROJECT_SOURCE_DIR}/Examples/CI/FO_Z_Variations/Base.yaml
            ${PROJECT_SOURCE_DIR}/Examples/CI/FO_Z_Variations/${order}.yaml
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
        )
        set_tests_properties(Main${test}partC${NN} PROPERTIES DEPENDS Main${test}partA)
        set_tests_properties(Main${test}partC${NN} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")
      endmacro()
      run_var(1 1.00 1.00 1   1)
      run_var(2 4.00 1.00 2   1)
      run_var(3 4.00 4.00 2   2)
      run_var(4 1.00 4.00 1   2)
      run_var(5 0.25 1.00 0.5 1)
      run_var(6 0.25 0.25 0.5 0.5)
      run_var(7 1.00 0.25 1   0.5)

      # finally, run python script to check for any diffs between the on-the-fly
      # and the explicit runs' analysis results
      add_test(
        NAME Main${test}partD
        COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/Examples/CI/FO_Z_Variations/diff_variations.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
        )
      set_tests_properties(Main${test}partD PROPERTIES
        ENVIRONMENT "${STANDARDENVIRONMENT}"
        DEPENDS "Main${test}partB;Main${test}partC1;Main${test}partC2;Main${test}partC3;Main${test}partC4;Main${test}partC5;Main${test}partC6;Main${test}partC7")
    endmacro()
    declare_variations_test(LO)
    declare_variations_test(NLO)
  endif()
endif()


if(SHERPA_ENABLE_CATCH2)
  # Add unit tests
  add_subdirectory(Unit)
endif()

if(SHERPA_ENABLE_UFO)
  # Add UFO tests
  add_subdirectory(End2End/UFO)
endif()

if(SHERPA_ENABLE_COVERAGE)
    find_program(LCOV lcov)
    find_program(GENHTML genhtml)

    if(LCOV AND GENHTML)
        get_filename_component(BUILD_DIR_NAME ${CMAKE_BINARY_DIR} NAME)
        get_filename_component(INSTALL_DIR_NAME ${CMAKE_INSTALL_PREFIX} NAME)
        set(COVERAGE_EXCL '${BUILD_DIR_NAME}/*' '${INSTALL_DIR_NAME}/*' '*/Tests/*' '*/yaml-cpp/*')
        set(LCOV_FLAGS --quiet --ignore-errors inconsistent,inconsistent)
        if(CMAKE_SYSTEM_NAME MATCHES Darwin)
          set(LCOV_FLAGS ${LCOV_FLAGS},usage,gcov,gcov)
        else()
        endif()
        add_custom_target(coverage
            COMMENT "SHERPA: Running coverage analysis with gcov/lcov..."
            COMMAND ${LCOV} --directory ${CMAKE_SOURCE_DIR} --no-external --capture --output-file coverage.info ${LCOV_FLAGS}
            COMMAND ${LCOV} --remove coverage.info --output-file coverage.info ${LCOV_FLAGS} ${COVERAGE_EXCL} || true
            COMMAND ${GENHTML} coverage.info --output-directory html ${LCOV_FLAGS}
            COMMAND ${CMAKE_COMMAND} -E echo "SHERPA: Coverage report generated in Tests/html/index.html"
        )
        add_custom_target(cov-reset
            COMMENT "SHERPA: Resetting coverage data..."
            COMMAND find ${CMAKE_BINARY_DIR} -name "*.gcda" -delete
        )
    else()
        message(WARNING "lcov or genhtml not found, coverage target will be unavailable.")
    endif()
endif()
