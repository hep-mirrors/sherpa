set(DireLorentz_SOURCES FFV_FF.C
                    FFV_FI.C
                    FFV_IF.C
                    FFV_II.C
                    VVV_FF.C
                    VVV_FI.C
                    VVV_IF.C
                    VVV_II.C
                    FFFF_FF.C
                    FFFF_FI.C
                    FFFF_IF.C
                    FFFF_II.C)

sherpa_create_git_info(DireLorentz ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DireLorentz SHARED ${DireLorentz_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DireLorentz DireLorentz)
add_dependencies(DireLorentz DireLorentz_Git_Info)
