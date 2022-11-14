##### Automatically or semiautomaticaly generated, do not edit.#########
#Loops/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsLoops.la)
 set(GITTAG  METOOLS/Loops)
set(Loops_SOURCES  Master_Tadpole.C Master_Bubble.C Master_Triangle.C Master_Box.C PV_Tadpoles.C PV_Bubbles.C PV_Triangles.C)
 #localincdir = $(pkg#includedir)/METOOLS/Loops)
set(Loops_HEADERS Divergence_Array.H Master_Integrals.H PV_Integrals.H)
creategitinfo(Loops ../Loops/)
list(TRANSFORM Loops_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Loops/")
list(APPEND Loops_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsLoops SHARED ${Loops_SOURCES} )
target_include_directories(METoolsLoops PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(METoolsLoops)
install(TARGETS METoolsLoops DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsLoops        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsLoops SOVERSION 0.0.0)
set_target_properties(METoolsLoops        PROPERTIES DEFINE_SYMBOL "")
