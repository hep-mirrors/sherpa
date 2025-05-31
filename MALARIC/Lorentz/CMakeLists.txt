set(MalaricLorentz_SOURCES Soft_FS.C
                    Soft_IS.C
                    FFV_FS.C
                    FFV_IS.C
                    VVV_FS.C
                    VVV_IS.C)

sherpa_create_git_info(MalaricLorentz ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(MalaricLorentz SHARED ${MalaricLorentz_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(MalaricLorentz MalaricLorentz)
