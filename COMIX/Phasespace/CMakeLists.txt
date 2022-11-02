
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND Phasespace_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ComixPhasespace SHARED 
                             ${Phasespace_SOURCES} 
)
target_include_directories(ComixPhasespace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ComixPhasespace DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ComixPhasespace        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ComixPhasespace SOVERSION 0.0.0)
set_target_properties(ComixPhasespace        PROPERTIES DEFINE_SYMBOL "")
