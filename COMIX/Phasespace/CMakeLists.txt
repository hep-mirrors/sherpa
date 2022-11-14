##### Automatically or semiautomaticaly generated, do not edit.#########
#Phasespace/Makefile.am
 
set(GITTAG  COMIX/Phasespace)
set(Phasespace_SOURCES  PS_Info.C PS_Vertex.C PS_Current.C PS_Generator.C PS_Channel.C)
set(Phasespace_HEADERS PS_Info.H PS_Vertex.H PS_Current.H PS_Generator.H PS_Channel.H)
sherpa_create_git_info(Phasespace ../Phasespace/)
list(TRANSFORM Phasespace_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Phasespace/")
list(APPEND Phasespace_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ComixPhasespace SHARED ${Phasespace_SOURCES} )
sherpa_handle_shared_library(ComixPhasespace ComixPhasespace)
