set(PhotonSplitter_SOURCES Kinematics_Base.C
                           Photon_Splitter.C
                           Splitting_Functions.C
                           Sudakov.C)

set(PhotonSplitter_HEADERS Kinematics_Base.H
                           Photon_Splitter.H
                           Splitting_Functions.H
                           Sudakov.H)


sherpa_create_git_info(PhotonSplitter ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PhotonsPhotonSplitter SHARED ${PhotonSplitter_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhotonsPhotonSplitter PhotonsPhotonSplitter)
add_dependencies(PhotonsPhotonSplitter PhotonSplitter_Git_Info)
