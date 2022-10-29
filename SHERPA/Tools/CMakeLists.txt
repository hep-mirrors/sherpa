
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaTools.la)
 set(GITTAG  SHERPA/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Output_Base.C Event_Reader_Base.C HepEvt_Interface.C HepMC2_Interface.C HepMC3_Interface.C Pythia_HepEvt_Translator.C Output_HepEvt.C Output_LHEF.C Userhook_Base.C Analysis_Interface.C)
set(Tools_HEADERS Definitions.H Output_Base.H Event_Reader_Base.H HepEvt_Interface.H HepMC2_Interface.H HepMC3_Interface.H Pythia_HepEvt_Translator.H Output_HepEvt.H Output_LHEF.H Userhook_Base.H Analysis_Interface.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaTools_static STATIC ${SherpaTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(SherpaTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaTools)
set_target_properties(SherpaTools_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaTools SHARED ${SherpaTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(SherpaTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaTools SOVERSION 0.0.0)
set_target_properties(SherpaTools        PROPERTIES DEFINE_SYMBOL "")
endif()
