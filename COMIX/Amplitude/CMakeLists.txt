
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Amplitude/Makefile.am
 
 #pkglib_LTLIBRARIES = libComixAmplitude.la)
 set(GITTAG  COMIX/Amplitude)
 #include ../../git.make)
set(Amplitude_SOURCES  Amplitude.C)
set(Amplitude_HEADERS Amplitude.H)
creategitinfo(Amplitude ../Amplitude/)
list(TRANSFORM Amplitude_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Amplitude/")
list(APPEND Amplitude_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ComixAmplitude SHARED 
                             ${Amplitude_SOURCES} 
)
target_include_directories(ComixAmplitude PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ComixAmplitude)
install(TARGETS ComixAmplitude DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ComixAmplitude        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ComixAmplitude SOVERSION 0.0.0)
set_target_properties(ComixAmplitude        PROPERTIES DEFINE_SYMBOL "")
