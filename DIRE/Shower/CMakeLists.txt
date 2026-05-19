set(DireShower_SOURCES Gauge.C
                   Lorentz.C
                   Kernel.C
                   Alpha_QCD.C
                   Lorentz_FF.C
                   Lorentz_FI.C
                   Lorentz_IF.C
                   Lorentz_II.C
                   Cluster_Definitions.C
                   Shower.C)

set(DireShower_HEADERS Gauge.H
                   Lorentz.H
                   Kernel.H
                   Alpha_QCD.H
                   Lorentz_FF.H
                   Lorentz_FI.H
                   Lorentz_IF.H
                   Lorentz_II.H
                   Cluster_Definitions.H
                   Shower.H)

sherpa_create_git_info(DireShower ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DireShower SHARED ${DireShower_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireShower DireShower)
add_dependencies(DireShower DireShower_Git_Info)
