set(YFSNLO_SOURCES NLO_Base.C
                 Real.C
                 Virtual.C)

set(YFSNLO_HEADERS  NLO_Base.H
                 Real.H
                 Virtual.H)


sherpa_create_git_info(YFSNLO)
 
add_library(YFSNLO SHARED ${YFSNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSNLO YFSNLO)
add_dependencies(YFSNLO YFSNLO_Git_Info)
