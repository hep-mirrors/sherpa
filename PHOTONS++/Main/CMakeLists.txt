set(PhotonsMain_SOURCES Photons.C
                 Define_Dipole.C
                 Dipole_Type.C)

set(PhotonsMain_HEADERS Photons.H
                 Dipole_Type.H
                 Define_Dipole.H)

sherpa_create_git_info(PhotonsMain)
 
add_library(PhotonsMain SHARED ${PhotonsMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhotonsMain PhotonsMain)
