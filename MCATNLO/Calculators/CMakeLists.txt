
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND Calculators_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(MCatNLOCalculators SHARED 
                             ${Calculators_SOURCES} 
)
target_include_directories(MCatNLOCalculators PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(MCatNLOCalculators)
install(TARGETS MCatNLOCalculators DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOCalculators        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOCalculators SOVERSION 0.0.0)
set_target_properties(MCatNLOCalculators        PROPERTIES DEFINE_SYMBOL "")
