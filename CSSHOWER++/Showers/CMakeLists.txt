set(CSSShowers_SOURCES Splitting_Function_Base.C
                    Splitting_Function_Group.C
                    SF_Coupling.C
                    SF_Lorentz.C
                    Sudakov.C
                    Shower.C
                    Kinematics_Base.C)

set(CSSShowers_HEADERS Splitting_Function_Base.H
                    Splitting_Function_Group.H
                    SF_Coupling.H
                    SF_Lorentz.H
                    SF_Key.H
                    Sudakov.H
                    Shower.H
                    Kinematics_Base.H)

sherpa_create_git_info(CSSShowers ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(CSShowers SHARED ${CSSShowers_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CSShowers CSShowers)
add_dependencies(CSShowers CSSShowers_Git_Info)
