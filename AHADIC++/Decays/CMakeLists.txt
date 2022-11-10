
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Decays/Makefile.am
 
 #pkglib_LTLIBRARIES = libAhadicDecays.la)
 set(GITTAG  AHADIC++/Decays)
 #include ../../git.make)
set(Decays_SOURCES  Cluster_Decayer.C Cluster_Splitter.C)
 #localincdir = $(pkg#includedir)/AHADIC++/Decays)
set(Decays_HEADERS Cluster_Decayer.H Cluster_Splitter.H)
creategitinfo(Decays ../Decays/)
list(TRANSFORM Decays_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Decays/")
list(APPEND Decays_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AhadicDecays SHARED 
                             ${Decays_SOURCES} 
)
target_include_directories(AhadicDecays PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(AhadicDecays)
install(TARGETS AhadicDecays DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicDecays        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AhadicDecays SOVERSION 0.0.0)
set_target_properties(AhadicDecays        PROPERTIES DEFINE_SYMBOL "")
