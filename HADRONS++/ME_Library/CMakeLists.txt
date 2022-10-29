
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#ME_Library/Makefile.am
 
 #pkglib_LTLIBRARIES = libHadronsMEs.la)
 set(GITTAG  HADRONS++/ME_Library)
 #include ../../git.make)
set(ME_Library_SOURCES  HD_ME_Base.C Generic.C Current_ME.C Eta_Decay_MEs.C Four_Body_MEs.C Three_Body_MEs.C Two_Body_MEs.C B_K_Semileptonic.C B_KStar_Semileptonic.C)
set(ME_Library_HEADERS HD_ME_Base.H Generic.H Current_ME.H PseudoScalar_Decay_MEs.H Eta_Decay_MEs.H Vector_Decay_MEs.H Phi_Decay_MEs.H Tensor_Decay_MEs.H Baryon_SimpleDecay_MEs.H HeavyBaryon_SimpleTwoBody_Decay_MEs.H
 Four_Body_MEs.H Three_Body_MEs.H Two_Body_MEs.H B_K_Semileptonic.H B_KStar_Semileptonic.H)
creategitinfo(ME_Library ../ME_Library/)
  list(TRANSFORM ME_Library_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../ME_Library/")
 
set(ME_Library_esources )
if (SHERPA_BUILD_STATIC)
add_library(HadronsMEs_static STATIC ${HadronsMEs_esources}
                             ${ME_Library_SOURCES} 
)
target_include_directories(HadronsMEs_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(HadronsMEs_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS HadronsMEs_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(HadronsMEs_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME HadronsMEs)
set_target_properties(HadronsMEs_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(HadronsMEs SHARED ${HadronsMEs_esources}
                             ${ME_Library_SOURCES} 
)
target_include_directories(HadronsMEs PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(HadronsMEs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS HadronsMEs DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(HadronsMEs        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME HadronsMEs SOVERSION 0.0.0)
set_target_properties(HadronsMEs        PROPERTIES DEFINE_SYMBOL "")
endif()
