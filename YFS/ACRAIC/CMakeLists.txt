set(ACRAIC_SOURCES Acraic_Amplitude.C
                 Acraic_Base.C
                 Acraic_Spinors.C
                 Acraic_Born.C
                 Acraic_Virtual.C
                 Acraic_Real.C
                 Acraic_Partitions.C
                 Acraic_Comix.C
                 Acraic_Checks.C
                 Acraic_Beta2.C
                 Acraic_LoopHelicity.C
                 Acraic_RealVirtualHelicity.C
                 Acraic_Resonance.C)

set(ACRAIC_HEADERS Acraic_Amplitude.H
                   Acraic_Base.H)


sherpa_create_git_info(ACRAIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(YFSACRAIC ${SHERPA_SUBLIB_TYPE} ${ACRAIC_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSACRAIC YFSACRAIC)
add_dependencies(YFSACRAIC ACRAIC_Git_Info)
