set(PS_Library_SOURCES HD_PS_Base.C
                       Three_Body_PSs.C
                       Two_Body_PSs.C
                       Four_Body_PSs.C
                       ResonanceFlavour.C)

set(PS_Library_HEADERS HD_PS_Base.H
                       Three_Body_PSs.H
                       Two_Body_PSs.H
                       Four_Body_PSs.H
                       ResonanceFlavour.H)

sherpa_create_git_info(PS_Library)
 
add_library(HadronsPSs SHARED ${PS_Library_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronsPSs HadronsPSs)
