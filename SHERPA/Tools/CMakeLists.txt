#Manually edited
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaTools.la)
 set(GITTAG  SHERPA/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Output_Base.C Event_Reader_Base.C HepEvt_Interface.C   Pythia_HepEvt_Translator.C Output_HepEvt.C Output_LHEF.C Userhook_Base.C Analysis_Interface.C)
if (SHERPA_ENABLE_HEPMC2)
  list(APPEND Tools_SOURCES HepMC2_Interface.C)
endif()
if (SHERPA_ENABLE_HEPMC3)
  list(APPEND Tools_SOURCES HepMC3_Interface.C)
endif()

set(Tools_HEADERS Definitions.H Output_Base.H Event_Reader_Base.H HepEvt_Interface.H HepMC2_Interface.H HepMC3_Interface.H Pythia_HepEvt_Translator.H Output_HepEvt.H Output_LHEF.H Userhook_Base.H Analysis_Interface.H)

set(foo_SOURCES  Output_Base.C Event_Reader_Base.C HepEvt_Interface.C HepMC2_Interface.C HepMC3_Interface.C Pythia_HepEvt_Translator.C Output_HepEvt.C Output_LHEF.C Userhook_Base.C Analysis_Interface.C)
set(foo_HEADERS Definitions.H Output_Base.H Event_Reader_Base.H HepEvt_Interface.H HepMC2_Interface.H HepMC3_Interface.H Pythia_HepEvt_Translator.H Output_HepEvt.H Output_LHEF.H Userhook_Base.H Analysis_Interface.H)

creategitinfo(foo ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaTools SHARED 
                             ${Tools_SOURCES} 
)
target_include_directories(SherpaTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SherpaTools)
install(TARGETS SherpaTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaTools SOVERSION 0.0.0)
set_target_properties(SherpaTools        PROPERTIES DEFINE_SYMBOL "")

if (SHERPA_ENABLE_HEPMC2)
  target_link_libraries(SherpaTools PRIVATE ${HEPMC2_LIBRARIES})
  target_include_directories(SherpaTools PRIVATE ${HEPMC2_INCLUDE_DIRS})
endif()
if (SHERPA_ENABLE_HEPMC3)
  target_link_libraries(SherpaTools PRIVATE ${HEPMC3_LIB})
  target_include_directories(SherpaTools PRIVATE ${HEPMC3_INCLUDE_DIR})
  if (SHERPA_ENABLE_HEPMC3_ROOT)
    target_link_libraries(SherpaTools PRIVATE ${HEPMC3_ROOTIO_LIB})
  endif()
endif()

if (SHERPA_ENABLE_HEPMC2)
target_link_libraries(SherpaTools PRIVATE ${HEPMC2_LIBRARIES})
target_include_directories(SherpaTools PRIVATE ${HEPMC2_INCLUDE_DIRS})
endif()
