
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND SMEHC_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSMEHC SHARED 
                             ${SMEHC_SOURCES} 
)
target_include_directories(SherpaSMEHC PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSMEHC DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSMEHC        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSMEHC SOVERSION 0.0.0)
set_target_properties(SherpaSMEHC        PROPERTIES DEFINE_SYMBOL "")
