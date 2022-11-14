##### Automatically or semiautomaticaly generated, do not edit.#########
#Shower/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireShower.la)
 set(GITTAG  DIRE/Shower)
set(Shower_SOURCES  Gauge.C Lorentz.C Kernel.C Alpha_QCD.C Lorentz_FF.C Lorentz_FI.C Lorentz_IF.C Lorentz_II.C Cluster_Definitions.C Shower.C)
set(Shower_HEADERS Gauge.H Lorentz.H Kernel.H Alpha_QCD.H Lorentz_FF.H Lorentz_FI.H Lorentz_IF.H Lorentz_II.H Cluster_Definitions.H Shower.H)
creategitinfo(Shower ../Shower/)
list(TRANSFORM Shower_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Shower/")
list(APPEND Shower_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DireShower SHARED ${Shower_SOURCES} )
target_include_directories(DireShower PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(DireShower)
install(TARGETS DireShower DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireShower        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireShower SOVERSION 0.0.0)
set_target_properties(DireShower        PROPERTIES DEFINE_SYMBOL "")
