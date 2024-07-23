set(NLO_SOURCES NLO_Base.C
                 Real.C
                 Virtual.C)

set(NLO_HEADERS  NLO_Base.H
                 Real.H
                 Virtual.H)


sherpa_create_git_info(NLO)
 
add_library(YFSNLO SHARED ${NLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSNLO YFSNLO)
