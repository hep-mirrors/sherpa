
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaMain.la)
 set(GITTAG  SHERPA/Main)
 #include ../../git.make)
set(Main_SOURCES  Sherpa.C Filter.C)
set(Main_HEADERS Sherpa.H Filter.H)
 #SYSLIBS = $(FCLIBS) -lm -ldl $(CONDITIONAL_THREADLIBS)
#libSherpaMain_la_LIBADD = @SHERPALIBS@ @SHRIMPSLIBS@ @HADRONSLIBS@ @RECONNECTIONSLIBS@ @AHADICLIBS@ @AMISICLIBS@ @PHOTONSLIBS@ @REMNANTSLIBS@ @MCATNLOLIBS@ @DIMLIBS@ @CSSLIBS@ @DIRELIBS@ @COMIXLIBS@ @AMEGICLIBS@ @EXTRAXSLIBS@ @EXTAMPLIBS@ @PHASICLIBS@ @METOOLSLIBS@ @MODELLIBS@ @PDFLIBS@ @BEAMLIBS@ @ATOOLSLIBS@ @CONDITIONAL_HEPMC3LIBS@ @CONDITIONAL_HEPMC2LIBS@ $(#SYSLIBS)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaMain_static STATIC ${SherpaMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(SherpaMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaMain)
set_target_properties(SherpaMain_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaMain SHARED ${SherpaMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(SherpaMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaMain SOVERSION 0.0.0)
set_target_properties(SherpaMain        PROPERTIES DEFINE_SYMBOL "")
endif()
