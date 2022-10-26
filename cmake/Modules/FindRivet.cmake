# - Try to find RIVET
# Defines:
#
#  RIVET_FOUND
#  RIVET_INCLUDE_DIR
#  RIVET_INCLUDE_DIRS (not cached)
#  RIVET_LIBRARY
#  RIVET_LIBRARIES (not cached)
#  RIVET_LIBRARY_DIRS (not cached)

if (RIVET_ROOT_DIR OR RIVET_DIR OR (DEFINED ENV{RIVET_ROOT_DIR}) OR (DEFINED ENV{RIVET_DIR}) )
  set(RIVET_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (RIVET_ROOT_DIR)
    list (APPEND RIVET_SEARCH_DIRS "${RIVET_ROOT_DIR}" )
  endif()
  if (RIVET_DIR)
    list (APPEND RIVET_SEARCH_DIRS "${RIVET_DIR}" )
  endif()
  if (DEFINED EVN{RIVET_ROOT_DIR})
    list (APPEND RIVET_SEARCH_DIRS "$ENV{RIVET_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{RIVET_DIR})
    list (APPEND RIVET_SEARCH_DIRS "ENV{RIVET_DIR}" )
  endif()
endif()

if (RIVET_SEARCH_DIRS)
  find_program(RIVET_EXE NAMES rivet PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES bin NO_DEFAULT_PATH )
  find_path(RIVET_DATA_PATH ATLAS_2012_I1094061.yoda PATH_SUFFIXES share/Rivet/)
  find_path(RIVET_INCLUDE_DIR Rivet/Rivet.hh PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_library(RIVET_LIBRARY NAMES Rivet PATHS ${RIVET_SEARCH_DIRS}  PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_program(RIVET_EXE NAMES rivet PATHS ${RIVET_SEARCH_DIRS} PATH_SUFFIXES bin)
  find_path(RIVET_DATA_PATH ATLAS_2012_I1094061.yoda PATH_SUFFIXES share/Rivet/)
  find_path(RIVET_INCLUDE_DIR Rivet/Rivet.hh PATH_SUFFIXES include)
  find_library(RIVET_LIBRARY NAMES Rivet PATHS_SUFFIXES lib lib64)
endif()
set(RIVET_VERSION 0.0.0)
if (RIVET_INCLUDE_DIR)
  if (EXISTS ${RIVET_INCLUDE_DIR}/Rivet/Config/RivetConfig.hh)
    file(STRINGS ${RIVET_INCLUDE_DIR}/Rivet/Config/RivetConfig.hh RIVET_VERSION_STRING_CONTENT REGEX "^#define[ ]+RIVET_VERSION[ ]+\"" )
    if (RIVET_VERSION_STRING_CONTENT)
      string(REGEX MATCH "[1234567890.]+[a-zA-Z]*" RIVET_VERSION ${RIVET_VERSION_STRING_CONTENT})
    endif()
    file(STRINGS ${RIVET_INCLUDE_DIR}/Rivet/Config/RivetConfig.hh RIVET_HEPMC_VERSION_STRING REGEX "^#define[ ]+RIVET_ENABLE_HEPMC_3[ ]+" )
    if (RIVET_HEPMC_VERSION_STRING)
      set(Rivet_HEPMC2_FOUND FALSE)
      set(Rivet_HEPMC3_FOUND TRUE)
    else()  
      set(Rivet_HEPMC2_FOUND TRUE)
      set(Rivet_HEPMC3_FOUND FALSE)
    endif()
  endif()
endif()


mark_as_advanced(RIVET_INCLUDE_DIR RIVET_LIBRARY RIVET_EXE)

# handle the QUIETLY and REQUIRED arguments and set RIVET_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Rivet HANDLE_COMPONENTS REQUIRED_VARS RIVET_INCLUDE_DIR RIVET_LIBRARY)

set(RIVET_LIBRARIES ${RIVET_LIBRARY})
get_filename_component(RIVET_LIBRARY_DIRS ${RIVET_LIBRARY} PATH)
get_filename_component(RIVET_ANALYSIS_PATH ${RIVET_LIBRARY} PATH)
set(RIVET_ANALYSIS_PATH ${RIVET_ANALYSIS_PATH}/Rivet)

set(RIVET_INCLUDE_DIRS ${RIVET_INCLUDE_DIR})

mark_as_advanced(RIVET_FOUND)

macro(RivetBuildPlugin a b)
  find_package(FastJet REQUIRED)
  string(REPLACE " " ";" sources_split ${b})
  add_library(RivetAnalysis${a} SHARED  ${sources_split} )
  set_target_properties(RivetAnalysis${a} PROPERTIES PREFIX "")
  target_include_directories(RivetAnalysis${a} PRIVATE ${RIVET_INCLUDE_DIRS} ${FASTJET_INCLUDE_DIR})
  target_link_libraries(RivetAnalysis${a} PRIVATE ${RIVET_LIBRARY} ${FASTJET_LIBRARIES})
  set_target_properties(RivetAnalysis${a} PROPERTIES 
                                               ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${a}/$<0:>
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${a}/$<0:>
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${a}/$<0:>)
  
  if (Rivet_HEPMC3_FOUND)
    target_compile_definitions(RivetAnalysis${a} PRIVATE -DENABLE_HEPMC_3=true)
  endif()
endmacro(RivetBuildPlugin a b)

macro(RivetRunPlugin a b c other)
  string(REPLACE " " "," analist ${b})
  string(REPLACE " " ";" optlist ${other})
  add_custom_target(RivetRunPlugin${a} 
                                        COMMAND  ${RIVET_EXE}   --analysis-path-append ${RIVET_ANALYSIS_PATH} -a ${analist} ${c} ${optlist}
                                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    )
endmacro(RivetRunPlugin a b c other)

macro(RivetMakeHTML a b c)
  string(REPLACE " " "," flist ${c})
  #string(REPLACE " " ";" optlist ${other})
  add_custom_target(RivetMakeHTML${a} 
                                        COMMAND  ${RIVET_EXE}-mkhtml ${b} ${flist}
                                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    )
endmacro(RivetMakeHTML a b c)



