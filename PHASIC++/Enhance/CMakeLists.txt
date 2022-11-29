##### Automatically or semiautomaticaly generated, do not edit.#########
#Enhance/Makefile.am
 
set(GITTAG  PHASIC++/Enhance)

set(Enhance_SOURCES Enhance_Observable_Base.C
                    Variable_Enhance_Observable.C
                    Fastjet_Enhance_Observable.C)

set(Enhance_HEADERS Enhance_Observable_Base.H)
sherpa_create_git_info(Enhance ../Enhance/)
list(TRANSFORM Enhance_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Enhance/")
list(APPEND Enhance_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PhasicEnhance SHARED ${Enhance_SOURCES} )
sherpa_handle_shared_library(PhasicEnhance PhasicEnhance)
