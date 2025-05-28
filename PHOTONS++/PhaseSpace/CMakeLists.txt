set(PhaseSpace_SOURCES Avarage_Photon_Number.C
                       Generate_Dipole_Photon_Angle.C
                       Generate_Multipole_Photon_Angle.C)

set(PhaseSpace_HEADERS Avarage_Photon_Number.H
                       Generate_Dipole_Photon_Angle.H
                       Generate_Multipole_Photon_Angle.H)

sherpa_create_git_info(PhaseSpace ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PhotonsPhaseSpace SHARED ${PhaseSpace_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhotonsPhaseSpace PhotonsPhaseSpace)
add_dependencies(PhotonsPhaseSpace PhaseSpace_Git_Info)
