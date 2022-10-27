
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Phasespace/Makefile.am
 
 #pkglib_LTLIBRARIES = libComixPhasespace.la)
 set(GITTAG  COMIX/Phasespace)
 #include ../../git.make)
set(Phasespace_SOURCES  PS_Info.C PS_Vertex.C PS_Current.C PS_Generator.C PS_Channel.C)
set(Phasespace_HEADERS PS_Info.H PS_Vertex.H PS_Current.H PS_Generator.H PS_Channel.H)
creategitinfo(Phasespace ../Phasespace/)
  list(TRANSFORM Phasespace_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Phasespace/")
 
set(Phasespace_esources )
if (SHERPA_BUILD_STATIC)
add_library(ComixPhasespace_static STATIC ${ComixPhasespace_esources}
                             ${Phasespace_SOURCES} 
)
target_include_directories(ComixPhasespace_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ComixPhasespace_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ComixPhasespace_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ComixPhasespace_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ComixPhasespace)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ComixPhasespace SHARED ${ComixPhasespace_esources}
                             ${Phasespace_SOURCES} 
)
target_include_directories(ComixPhasespace PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ComixPhasespace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ComixPhasespace DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ComixPhasespace        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ComixPhasespace SOVERSION 0.0.0)
endif()
