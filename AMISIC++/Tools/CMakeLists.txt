##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libAmisicTools.la)
 set(GITTAG  AMISIC++/Tools)
set(Tools_SOURCES  MI_Parameters.C MPI_Scale_Setter.C MPI_KFactor_Setter.C Hadronic_XSec_Calculator.C Impact_Parameter.C Matter_Overlap.C Interaction_Probability.C Over_Estimator.C)
set(Tools_HEADERS MI_Parameters.H MPI_Scale_Setter.H MPI_KFactor_Setter.H Hadronic_XSec_Calculator.H Impact_Parameter.H Matter_Overlap.H Interaction_Probability.H Over_Estimator.H)
creategitinfo(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AmisicTools SHARED ${Tools_SOURCES} )
target_include_directories(AmisicTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(AmisicTools)
install(TARGETS AmisicTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AmisicTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AmisicTools SOVERSION 0.0.0)
set_target_properties(AmisicTools        PROPERTIES DEFINE_SYMBOL "")
