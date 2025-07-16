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
                WORKING_DIRECTORY ${WORKING_DIR}
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
string(REPLACE "," ";" ALL_FILES "${ALL_FILES_ARG}")
list (FILTER ALL_FILES EXCLUDE REGEX ".*CXXFLAGS.*")
set(content "")
foreach (f IN LISTS ALL_FILES )
    file(READ "${WORKING_DIR}/${f}"  temp)
  string(APPEND content "${temp}")
endforeach()
string(REPLACE "${PROJECT_SOURCE_DIR}/" "" GITTAG "${WORKING_DIR}")
string(MD5 MDF "${content}")
set(newgitinfo "#include \"ATOOLS/Org/Git_Info.H\"\nstatic ATOOLS::Git_Info initializer\n(\"${GITTAG}\",\"${GITURL}\",\"${GITREV}${GITREVSUFFIX}\",\"${MDF}\");\n")
set(oldgitinfo "")
if (EXISTS ${BINARYDIR}/Git_Info.C)
  file(READ ${BINARYDIR}/Git_Info.C oldgitinfo)
endif()
if ( NOT "${newgitinfo}" STREQUAL "${oldgitinfo}")
  file( WRITE ${BINARYDIR}/Git_Info.C "${newgitinfo}")
endif()
list(TRANSFORM ALL_FILES PREPEND "${WORKING_DIR}/")
set_property(GLOBAL APPEND PROPERTY ALL_SOURCE_FILES ${ALL_FILES})
