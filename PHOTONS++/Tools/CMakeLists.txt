##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  PHOTONS++/Tools)

set(Tools_SOURCES Dipole_FF.C
                  Dipole_FI.C
                  Dress_Blob_Base.C
                  Generate_One_Photon.C
                  Weight_Dipole.C
                  Weight_Higher_Order_Corrections.C
                  Weight_Jacobian.C
                  Weight_YFS.C
                  YFS_Form_Factor.C)

set(Tools_HEADERS Dipole_FF.H
                  Dipole_FI.H
                  Dress_Blob_Base.H
                  Generate_One_Photon.H
                  Weight_Base.H
                  Weight_Dipole.H
                  Weight_Higher_Order_Corrections.H
                  Weight_Jacobian.H
                  Weight_YFS.H
                  YFS_Form_Factor.H)

sherpa_create_git_info(Tools)
 
add_library(PhotonsTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PhotonsTools PhotonsTools)
