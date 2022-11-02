
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Shower/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMShower.la)
 set(GITTAG  DIM/Shower)
 #include ../../git.make)
set(Shower_SOURCES  Gauge.C Lorentz.C Kernel.C Alpha_QCD.C Lorentz_FF.C Lorentz_FI.C Lorentz_IF.C Lorentz_II.C Shower.C)
set(Shower_HEADERS Gauge.H Lorentz.H Kernel.H Alpha_QCD.H Lorentz_FF.H Lorentz_FI.H Lorentz_IF.H Lorentz_II.H Shower.H)
creategitinfo(Shower ../Shower/)
list(TRANSFORM Shower_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Shower/")
list(APPEND Shower_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DIMShower SHARED 
                             ${Shower_SOURCES} 
)
target_include_directories(DIMShower PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DIMShower DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMShower        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMShower SOVERSION 0.0.0)
set_target_properties(DIMShower        PROPERTIES DEFINE_SYMBOL "")
