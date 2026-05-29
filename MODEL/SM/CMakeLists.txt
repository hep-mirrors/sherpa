set(SM_SOURCES Model.C
               LowEnergy_Model.C 
               Baryon_FormFactors.C 
               Meson_FormFactors.C )

set(SM_HEADERS Model.H
               LowEnergy_Model.H )

sherpa_create_git_info(SM ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaSM SHARED ${SM_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaSM SherpaSM)
add_dependencies(SherpaSM SM_Git_Info)
