find_program(GIT NAMES git)
if (GIT AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(COMMAND ${GIT} rev-parse --abbrev-ref HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITURL
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GIT} rev-parse  HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITREV
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  message(STATUS "SHERPA: GIT IS NOT AVAILABLE!")
  set(GITURL "unknownurl")
  set(GITREV "unknownrevision")
endif()
if (GIT AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(COMMAND ${GIT} status -s --untracked-files=no .
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE GITCLEAN
                OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(GITCLEAN "")
endif()
if ("${GITCLEAN}" STREQUAL "")
  set(GITREVSUFFIX "")
else()
  set(GITREVSUFFIX "-dirty")
endif()
set(allsources )
separate_arguments(allsources UNIX_COMMAND "${ALL_SOURCES} ${ALL_HEADERS}")
list (FILTER allsources EXCLUDE REGEX ".*CXXFLAGS.*")
set(content "")
foreach (f IN LISTS allsources )
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/${f}"  temp)
  string(APPEND content "${temp}")
endforeach()
string(REPLACE "${PROJECT_SOURCE_DIR}/" "" GITTAG "${CMAKE_CURRENT_SOURCE_DIR}")
string(MD5 MDF "${content}")
set(newgitinfo "#include \"ATOOLS/Org/Git_Info.H\"\nstatic ATOOLS::Git_Info initializer\n(\"${GITTAG}\",\"${GITURL}\",\"${GITREV}${GITREVSUFFIX}\",\"${MDF}\");\n")
set(oldgitinfo "")
if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C oldgitinfo)
endif()
if ( NOT "${newgitinfo}" STREQUAL "${oldgitinfo}")
  file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C "${newgitinfo}")
endif()
list(TRANSFORM allsources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
set_property(GLOBAL APPEND PROPERTY ALL_SOURCE_FILES ${allsources})
