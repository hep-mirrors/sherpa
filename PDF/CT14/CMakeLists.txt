##### Automatically or semiautomaticaly generated, do not edit.#########
#CT14/Makefile.am
 
 #pkglib_LTLIBRARIES = libCT14Sherpa.la)
 set(GITTAG  PDF/CT14)
set(CT14_SOURCES  CT14Pdf.f CT14_Fortran_Interface.C)
#nobase_#dist_pkgdata_DATA = CT14Grid/ct14n.00.pds CT14Grid/ct14nn.00.pds
creategitinfo(CT14 ../CT14/)
list(TRANSFORM CT14_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../CT14/")
list(APPEND CT14_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(CT14Sherpa SHARED ${CT14_SOURCES} )
target_include_directories(CT14Sherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(CT14Sherpa)
install(TARGETS CT14Sherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CT14Sherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CT14Sherpa SOVERSION 0.0.0)
set_target_properties(CT14Sherpa        PROPERTIES DEFINE_SYMBOL "")
