
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libAhadicTools.la)
 set(GITTAG  AHADIC++/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Proto_Particle.C Cluster.C Constituents.C Wave_Function.C Multiplet_Constructor.C Single_Transitions.C Double_Transitions.C Hadronisation_Parameters.C Soft_Cluster_Handler.C Singlet_Tools.C
 Splitter_Base.C Flavour_Selector.C KT_Selector.C Z_Selector.C)
 #localincdir = $(pkg#includedir)/AHADIC++/Tools)
set(Tools_HEADERS Proto_Particle.H Cluster.H Constituents.H Wave_Function.H Hadron_Multiplet.H Multiplet_Constructor.H Single_Transitions.H Double_Transitions.H Hadronisation_Parameters.H Soft_Cluster_Handler.H
 Singlet_Tools.H Splitter_Base.H Flavour_Selector.H KT_Selector.H Z_Selector.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(AhadicTools_static STATIC ${AhadicTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(AhadicTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME AhadicTools)
set_target_properties(AhadicTools_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(AhadicTools SHARED ${AhadicTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(AhadicTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AhadicTools SOVERSION 0.0.0)
set_target_properties(AhadicTools        PROPERTIES DEFINE_SYMBOL "")
endif()
