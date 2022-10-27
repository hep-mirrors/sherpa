
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
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
 
set(Amplitude_esources )
if (SHERPA_BUILD_STATIC)
add_library(ComixAmplitude_static STATIC ${ComixAmplitude_esources}
                             ${Amplitude_SOURCES} 
)
target_include_directories(ComixAmplitude_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ComixAmplitude_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ComixAmplitude_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ComixAmplitude_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ComixAmplitude)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ComixAmplitude SHARED ${ComixAmplitude_esources}
                             ${Amplitude_SOURCES} 
)
target_include_directories(ComixAmplitude PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ComixAmplitude PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ComixAmplitude DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ComixAmplitude        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ComixAmplitude SOVERSION 0.0.0)
endif()
