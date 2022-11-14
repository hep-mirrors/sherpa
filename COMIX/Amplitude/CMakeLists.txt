##### Automatically or semiautomaticaly generated, do not edit.#########
#Amplitude/Makefile.am
 
 #pkglib_LTLIBRARIES = libComixAmplitude.la)
 set(GITTAG  COMIX/Amplitude)
set(Amplitude_SOURCES  Amplitude.C)
set(Amplitude_HEADERS Amplitude.H)
sherpa_create_git_info(Amplitude ../Amplitude/)
list(TRANSFORM Amplitude_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Amplitude/")
list(APPEND Amplitude_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ComixAmplitude SHARED ${Amplitude_SOURCES} )
sherpa_handle_shared_library(ComixAmplitude ComixAmplitude)
