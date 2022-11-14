##### Automatically or semiautomaticaly generated, do not edit.#########
#Calculators/Makefile.am
 
 #pkglib_LTLIBRARIES = libCSCalculators.la)
 set(GITTAG  CSSHOWER++/Calculators)
set(Calculators_SOURCES  CF_QCD.C CF_QED.C CF_EW_FFV.C CF_Triangle.C CF_HV.C LF_SSV.C LF_FFV.C LF_VVV.C LF_Triangle.C)
creategitinfo(Calculators ../Calculators/)
list(TRANSFORM Calculators_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Calculators/")
list(APPEND Calculators_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(CSCalculators SHARED ${Calculators_SOURCES} )
target_include_directories(CSCalculators PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(CSCalculators)
install(TARGETS CSCalculators DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSCalculators        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CSCalculators SOVERSION 0.0.0)
set_target_properties(CSCalculators        PROPERTIES DEFINE_SYMBOL "")
