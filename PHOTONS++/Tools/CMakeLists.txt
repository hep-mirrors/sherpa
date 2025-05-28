set(PhotonsTools_SOURCES Dipole_FF.C
                  Dipole_FI.C
                  Dress_Blob_Base.C
                  Generate_One_Photon.C
                  Weight_Dipole.C
                  Weight_Higher_Order_Corrections.C
                  Weight_Jacobian.C
                  Weight_YFS.C
                  YFS_Form_Factor.C)

set(PhotonsTools_HEADERS Dipole_FF.H
                  Dipole_FI.H
                  Dress_Blob_Base.H
                  Generate_One_Photon.H
                  Weight_Base.H
                  Weight_Dipole.H
                  Weight_Higher_Order_Corrections.H
                  Weight_Jacobian.H
                  Weight_YFS.H
                  YFS_Form_Factor.H)

sherpa_create_git_info(PhotonsTools ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PhotonsTools SHARED ${PhotonsTools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhotonsTools PhotonsTools)
add_dependencies(PhotonsTools PhotonsTools_Git_Info)
