set(Lorentz_SOURCES FFV_FF.C
                    FFV_FI.C
                    FFV_IF.C
                    FFV_II.C
                    VVV_FF.C
                    VVV_FI.C
                    VVV_IF.C
                    VVV_II.C)

sherpa_create_git_info(Lorentz)
 
add_library(DIMLorentz SHARED ${Lorentz_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DIMLorentz DIMLorentz)
