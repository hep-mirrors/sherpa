set(DIMShower_SOURCES Gauge.C
                   Lorentz.C
                   Kernel.C
                   Alpha_QCD.C
                   Lorentz_FF.C
                   Lorentz_FI.C
                   Lorentz_IF.C
                   Lorentz_II.C
                   Shower.C)

set(DIMShower_HEADERS Gauge.H
                   Lorentz.H
                   Kernel.H
                   Alpha_QCD.H
                   Lorentz_FF.H
                   Lorentz_FI.H
                   Lorentz_IF.H
                   Lorentz_II.H
                   Shower.H)

sherpa_create_git_info(DIMShower ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DIMShower SHARED ${DIMShower_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DIMShower DIMShower)
add_dependencies(DIMShower DIMShower_Git_Info)
