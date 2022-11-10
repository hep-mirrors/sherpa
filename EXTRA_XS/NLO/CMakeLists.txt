
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#NLO/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXSNLO.la)
 set(GITTAG  EXTRA_XS/NLO)
 #include ../../git.make)
set(NLO_SOURCES  Logarithms.C Dummy_QCD_Virtual.C DY_QCD_Virtual.C DY_Jet_QCD_Virtual.C Higgs_QCD_Virtual.C Higgs_Jet_QCD_Virtual.C ee3jet_CSS_approx.C ppll1jet_CSS_approx.C epe2jet_CSS_approx.C)
set(NLO_HEADERS Logarithms.H)
creategitinfo(NLO ../NLO/)
list(TRANSFORM NLO_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../NLO/")
list(APPEND NLO_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXSNLO SHARED 
                             ${NLO_SOURCES} 
)
target_include_directories(ExtraXSNLO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ExtraXSNLO)
install(TARGETS ExtraXSNLO DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXSNLO        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXSNLO SOVERSION 0.0.0)
set_target_properties(ExtraXSNLO        PROPERTIES DEFINE_SYMBOL "")
