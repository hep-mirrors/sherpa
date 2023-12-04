set(Shower_SOURCES Gauge.C
                   Lorentz.C
                   Kernel.C
                   Alpha_QCD.C
                   Lorentz_FS.C
                   Lorentz_IS.C
                   Cluster_Definitions.C
                   Shower.C)

set(Shower_HEADERS Gauge.H
                   Lorentz.H
                   Kernel.H
                   Alpha_QCD.H
                   Lorentz_FS.H
                   Lorentz_IS.H
                   Cluster_Definitions.H
                   Shower.H)

sherpa_create_git_info(Shower)
 
add_library(AlaricShower SHARED ${Shower_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AlaricShower AlaricShower)
