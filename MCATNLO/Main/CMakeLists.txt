
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOMain.la)
 set(GITTAG  MCATNLO/Main)
 #include ../../git.make)
set(Main_SOURCES  CS_Cluster_Definitions.C CS_Gamma.C CS_MCatNLO.C)
set(Main_HEADERS CS_Cluster_Definitions.H CS_Gamma.H CS_MCatNLO.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(MCatNLOMain SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(MCatNLOMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(MCatNLOMain)
install(TARGETS MCatNLOMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOMain SOVERSION 0.0.0)
set_target_properties(MCatNLOMain        PROPERTIES DEFINE_SYMBOL "")
