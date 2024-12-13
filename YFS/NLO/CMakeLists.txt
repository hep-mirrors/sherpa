set(YFSNLO_SOURCES NLO_Base.C
                 Real.C
                 RealReal.C
                 RealVirtual.C
                 NLO_Amp.C
                 Virtual.C)

set(YFSNLO_HEADERS  NLO_Base.H
                 Real.H
                 RealReal.H
                 RealVirtual.H
                 NLO_Amp.H
                 Virtual.H)


sherpa_create_git_info(YFSNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(YFSNLO SHARED ${YFSNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSNLO YFSNLO)
add_dependencies(YFSNLO YFSNLO_Git_Info)
