find_package(Catch2 3 REQUIRED)

# find all .C files in the unit tests directory (and subdirectories)
file(GLOB_RECURSE UNIT_TEST_SRC "${CMAKE_SOURCE_DIR}/Tests/Unit/test_*.C")
add_executable(sherpa_tests "${UNIT_TEST_SRC}")

# add the additional include directories to the sherpa_tests target
target_include_directories(sherpa_tests PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
)

if(NOT (CMAKE_SYSTEM_NAME MATCHES Darwin))
  target_link_options(sherpa_tests PRIVATE "LINKER:-no-as-needed")
endif()

# for Catch2 v2 change Catch2::Catch2WithMain to Catch2::Catch2
target_link_libraries(sherpa_tests PRIVATE Catch2::Catch2WithMain SHERPAAll)

#add_test(unit_tests sherpa_tests)
include(CTest)
include(Catch)
catch_discover_tests(sherpa_tests)

