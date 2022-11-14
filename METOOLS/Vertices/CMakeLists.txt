##### Automatically or semiautomaticaly generated, do not edit.#########
#Vertices/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsVertices.la)
 set(GITTAG  METOOLS/Vertices)
set(Vertices_SOURCES  FFV_DLC.C VVV_DLC.C SSS_LC.C SSV_LC.C SFF_LC.C SVV_LC.C FFV_LC.C VVV_LC.C VVP_LC.C SSSS_LC.C SSVV_LC.C VVVV_LC.C)
creategitinfo(Vertices ../Vertices/)
list(TRANSFORM Vertices_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Vertices/")
list(APPEND Vertices_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsVertices SHARED ${Vertices_SOURCES} )
target_include_directories(METoolsVertices PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(METoolsVertices)
install(TARGETS METoolsVertices DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsVertices        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsVertices SOVERSION 0.0.0)
set_target_properties(METoolsVertices        PROPERTIES DEFINE_SYMBOL "")
