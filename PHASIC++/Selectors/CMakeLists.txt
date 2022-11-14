##### Automatically or semiautomaticaly generated, do not edit.#########
#Selectors/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicSelectors.la)
 set(GITTAG  PHASIC++/Selectors)
set(Selectors_SOURCES  Cut_Data.C Selector.C Variable_Selector.C Combined_Selector.C FastJet_Criterion.C Jet_Finder.C Jet_Selector.C Isolation_Selector.C NJet_Finder.C Standard_Selector.C DIS_Selector.C
 Decay_Selector.C MinSelector.C DressedParticleSelector.C Fastjet_Selector_Base.C Fastjet_Finder.C Fastjet_Selector.C Fastjet_Veto.C)
set(Selectors_HEADERS Cut_Data.H Selector.H Combined_Selector.H Jet_Finder.H Fastjet_Selector_Base.H)
sherpa_create_git_info(Selectors ../Selectors/)
list(TRANSFORM Selectors_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Selectors/")
list(APPEND Selectors_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PhasicSelectors SHARED ${Selectors_SOURCES} )
sherpa_handle_shared_library(PhasicSelectors PhasicSelectors)
