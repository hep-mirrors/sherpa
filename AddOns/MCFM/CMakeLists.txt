add_library(SherpaMCFM SHARED MCFM_Interface.C)
install(TARGETS SherpaMCFM DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaMCFM       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaMCFM SOVERSION 0.0.0)
target_include_directories(SherpaMCFM  PRIVATE ${MCFM_INCLUDE_DIR})
target_link_libraries(SherpaMCFM  PRIVATE ${MCFM_LIBRARY})
set_target_properties(SherpaMCFM       PROPERTIES DEFINE_SYMBOL "")
target_compile_options(SherpaMCFM       PRIVATE -fopenmp)
#-fopenmp ?
