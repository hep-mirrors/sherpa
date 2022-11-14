##### Automatically or semiautomaticaly generated, do not edit.#########
#Shower/Makefile.am
 
set(GITTAG  DIM/Shower)
set(Shower_SOURCES  Gauge.C Lorentz.C Kernel.C Alpha_QCD.C Lorentz_FF.C Lorentz_FI.C Lorentz_IF.C Lorentz_II.C Shower.C)
set(Shower_HEADERS Gauge.H Lorentz.H Kernel.H Alpha_QCD.H Lorentz_FF.H Lorentz_FI.H Lorentz_IF.H Lorentz_II.H Shower.H)
sherpa_create_git_info(Shower ../Shower/)
list(TRANSFORM Shower_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Shower/")
list(APPEND Shower_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(DIMShower SHARED ${Shower_SOURCES} )
sherpa_handle_shared_library(DIMShower DIMShower)
