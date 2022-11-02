
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Selectors/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicSelectors.la)
 set(GITTAG  PHASIC++/Selectors)
 #include ../../git.make)
set(Selectors_SOURCES  Cut_Data.C Selector.C Variable_Selector.C Combined_Selector.C FastJet_Criterion.C Jet_Finder.C Jet_Selector.C Isolation_Selector.C NJet_Finder.C Standard_Selector.C DIS_Selector.C
 Decay_Selector.C MinSelector.C DressedParticleSelector.C Fastjet_Selector_Base.C Fastjet_Finder.C Fastjet_Selector.C Fastjet_Veto.C)
set(Selectors_HEADERS Cut_Data.H Selector.H Combined_Selector.H Jet_Finder.H Fastjet_Selector_Base.H)
creategitinfo(Selectors ../Selectors/)
list(TRANSFORM Selectors_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Selectors/")
list(APPEND Selectors_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PhasicSelectors SHARED 
                             ${Selectors_SOURCES} 
)
target_include_directories(PhasicSelectors PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicSelectors DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicSelectors        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhasicSelectors SOVERSION 0.0.0)
set_target_properties(PhasicSelectors        PROPERTIES DEFINE_SYMBOL "")
