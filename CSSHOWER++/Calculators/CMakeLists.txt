
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Calculators/Makefile.am
 
 #pkglib_LTLIBRARIES = libCSCalculators.la)
 set(GITTAG  CSSHOWER++/Calculators)
 #include ../../git.make)
set(Calculators_SOURCES  CF_QCD.C CF_QED.C CF_EW_FFV.C CF_Triangle.C CF_HV.C LF_SSV.C LF_FFV.C LF_VVV.C LF_Triangle.C)
creategitinfo(Calculators ../Calculators/)
  list(TRANSFORM Calculators_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Calculators/")
 
set(Calculators_esources )
if (SHERPA_BUILD_STATIC)
add_library(CSCalculators_static STATIC ${CSCalculators_esources}
                             ${Calculators_SOURCES} 
)
target_include_directories(CSCalculators_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSCalculators_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSCalculators_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSCalculators_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME CSCalculators)
endif()
if (SHERPA_BUILD_SHARED)
add_library(CSCalculators SHARED ${CSCalculators_esources}
                             ${Calculators_SOURCES} 
)
target_include_directories(CSCalculators PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSCalculators PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSCalculators DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSCalculators        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CSCalculators SOVERSION 0.0.0)
endif()
