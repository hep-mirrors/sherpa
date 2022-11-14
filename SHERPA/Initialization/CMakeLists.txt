##### Automatically or semiautomaticaly generated, do not edit.#########
#Initialization/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaInitialization.la)
 set(GITTAG  SHERPA/Initialization)
set(Initialization_SOURCES  Initialization_Handler.C)
set(Initialization_HEADERS Initialization_Handler.H)
creategitinfo(Initialization ../Initialization/)
list(TRANSFORM Initialization_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Initialization/")
list(APPEND Initialization_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaInitialization SHARED ${Initialization_SOURCES} )
target_include_directories(SherpaInitialization PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SherpaInitialization)
install(TARGETS SherpaInitialization DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaInitialization        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaInitialization SOVERSION 0.0.0)
set_target_properties(SherpaInitialization        PROPERTIES DEFINE_SYMBOL "")
