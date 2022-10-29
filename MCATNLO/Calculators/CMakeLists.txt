
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Calculators/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOCalculators.la)
 set(GITTAG  MCATNLO/Calculators)
 #include ../../git.make)
set(Calculators_SOURCES  CF_QCD.C LF_FFV.C LF_VVV.C)
creategitinfo(Calculators ../Calculators/)
  list(TRANSFORM Calculators_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Calculators/")
 
set(Calculators_esources )
if (SHERPA_BUILD_STATIC)
add_library(MCatNLOCalculators_static STATIC ${MCatNLOCalculators_esources}
                             ${Calculators_SOURCES} 
)
target_include_directories(MCatNLOCalculators_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOCalculators_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOCalculators_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOCalculators_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME MCatNLOCalculators)
set_target_properties(MCatNLOCalculators_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(MCatNLOCalculators SHARED ${MCatNLOCalculators_esources}
                             ${Calculators_SOURCES} 
)
target_include_directories(MCatNLOCalculators PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOCalculators PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOCalculators DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOCalculators        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOCalculators SOVERSION 0.0.0)
set_target_properties(MCatNLOCalculators        PROPERTIES DEFINE_SYMBOL "")
endif()
