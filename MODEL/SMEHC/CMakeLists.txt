
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#SMEHC/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaSMEHC.la)
 set(GITTAG  MODEL/SMEHC)
 #include ../../git.make)
set(SMEHC_SOURCES  Effective_Higgs_Coupling.C Model.C SVV_LC.C SVVV_LC.C SVVP_LC.C SVVVV_LC.C)
set(SMEHC_HEADERS Effective_Higgs_Coupling.H)
creategitinfo(SMEHC ../SMEHC/)
  list(TRANSFORM SMEHC_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SMEHC/")
 
set(SMEHC_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaSMEHC_static STATIC ${SherpaSMEHC_esources}
                             ${SMEHC_SOURCES} 
)
target_include_directories(SherpaSMEHC_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSMEHC_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSMEHC_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSMEHC_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaSMEHC)
set_target_properties(SherpaSMEHC_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaSMEHC SHARED ${SherpaSMEHC_esources}
                             ${SMEHC_SOURCES} 
)
target_include_directories(SherpaSMEHC PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSMEHC PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSMEHC DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSMEHC        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSMEHC SOVERSION 0.0.0)
set_target_properties(SherpaSMEHC        PROPERTIES DEFINE_SYMBOL "")
endif()
