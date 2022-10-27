
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
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
 
set(Decays_esources )
if (SHERPA_BUILD_STATIC)
add_library(AhadicDecays_static STATIC ${AhadicDecays_esources}
                             ${Decays_SOURCES} 
)
target_include_directories(AhadicDecays_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicDecays_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicDecays_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicDecays_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME AhadicDecays)
endif()
if (SHERPA_BUILD_SHARED)
add_library(AhadicDecays SHARED ${AhadicDecays_esources}
                             ${Decays_SOURCES} 
)
target_include_directories(AhadicDecays PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicDecays PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicDecays DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicDecays        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AhadicDecays SOVERSION 0.0.0)
endif()
