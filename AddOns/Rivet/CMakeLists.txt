#Manually edited
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Rivet/Makefile.am
 
 set(GITTAG  AddOns/Rivet)
 #include ../../git.make)
set(Rivet_SOURCES  Rivet_Interface.C)
 if (SHERPA_ENABLE_RIVET)
 #pkglib_LTLIBRARIES = libSherpaRivetAnalysis.la)
set(Rivet_SOURCES ${Rivet_SOURCES})
 else()
 #RIVET_#EXTRA_DIST = $(RIVET_SOURCES)
 endif()
 #EXTRA_DIST = $(RIVET_#EXTRA_DIST)
creategitinfo(Rivet ../Rivet/)
list(TRANSFORM Rivet_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Rivet/")
list(APPEND Rivet_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaRivetAnalysis SHARED 
                             ${Rivet_SOURCES} 
)
target_include_directories(SherpaRivetAnalysis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaRivetAnalysis DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaRivetAnalysis        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaRivetAnalysis SOVERSION 0.0.0)
set_target_properties(SherpaRivetAnalysis        PROPERTIES DEFINE_SYMBOL "")
target_link_libraries(SherpaRivetAnalysis PRIVATE ${RIVET_LIBRARIES})
target_include_directories(SherpaRivetAnalysis PRIVATE ${RIVET_INCLUDE_DIR})
