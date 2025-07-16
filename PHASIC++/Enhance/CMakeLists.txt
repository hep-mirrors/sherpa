set(Enhance_SOURCES Enhance_Observable_Base.C
                    Variable_Enhance_Observable.C
                    Fastjet_Enhance_Observable.C)

set(Enhance_HEADERS Enhance_Observable_Base.H)

sherpa_create_git_info(Enhance ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PhasicEnhance SHARED ${Enhance_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhasicEnhance PhasicEnhance)
add_dependencies(PhasicEnhance Enhance_Git_Info)
