set(Shower_SOURCES Gauge.C
                   Lorentz.C
                   Kernel.C
                   Alpha_QCD.C
                   Lorentz_FF.C
                   Lorentz_FI.C
                   Lorentz_IF.C
                   Lorentz_II.C
                   Cluster_Definitions.C
                   Shower.C)

set(Shower_HEADERS Gauge.H
                   Lorentz.H
                   Kernel.H
                   Alpha_QCD.H
                   Lorentz_FF.H
                   Lorentz_FI.H
                   Lorentz_IF.H
                   Lorentz_II.H
                   Cluster_Definitions.H
                   Shower.H)

sherpa_create_git_info(Shower)
 
add_library(DireShower SHARED ${Shower_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireShower DireShower)
