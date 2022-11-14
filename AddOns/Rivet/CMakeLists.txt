#Manually edited
####  Automatically or semiautomaticaly generated, do not edit.
#Rivet/Makefile.am
 
 set(GITTAG  AddOns/Rivet)
set(Rivet_SOURCES  Rivet_Interface.C)
 if (SHERPA_ENABLE_RIVET)
set(Rivet_SOURCES ${Rivet_SOURCES})
 else()
 #RIVET_#EXTRA_DIST = $(RIVET_SOURCES)
 endif()
 #EXTRA_DIST = $(RIVET_#EXTRA_DIST)
sherpa_create_git_info(Rivet ../Rivet/)
list(TRANSFORM Rivet_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Rivet/")
list(APPEND Rivet_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaRivetAnalysis SHARED  ${Rivet_SOURCES} )
sherpa_handle_shared_library(SherpaRivetAnalysis SherpaRivetAnalysis)
target_link_libraries(SherpaRivetAnalysis PRIVATE ${RIVET_LIBRARIES})
target_include_directories(SherpaRivetAnalysis PRIVATE ${RIVET_INCLUDE_DIR} ${RIVET_CONFIG_CPPFLAGS_DIRS})
if (Rivet_HEPMC2_FOUND)
  target_include_directories(SherpaRivetAnalysis PRIVATE ${HEPMC2_INCLUDE_DIRS})
endif()
if (Rivet_HEPMC3_FOUND)
  target_include_directories(SherpaRivetAnalysis PRIVATE ${HEPMC3_INCLUDE_DIR})
endif()
