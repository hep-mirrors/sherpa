set(DipoleSubtraction_SOURCES DipoleSplitting_Base.C
                              FF_DipoleSplitting.C
                              FI_DipoleSplitting.C
                              IF_DipoleSplitting.C
                              II_DipoleSplitting.C
                              Single_LOProcess.C
                              Single_LOProcess_MHV.C
                              Single_LOProcess_External.C
                              Single_DipoleTerm.C
                              Single_OSTerm.C
                              Single_Virtual_Correction.C
                              Single_Real_Correction.C)

set(DipoleSubtraction_HEADERS DipoleSplitting_Base.H
                              FF_DipoleSplitting.H
                              FI_DipoleSplitting.H
                              IF_DipoleSplitting.H
                              II_DipoleSplitting.H
                              Single_LOProcess.H
                              Single_LOProcess_MHV.H
                              Single_LOProcess_External.H
                              Single_DipoleTerm.H
                              Single_OSTerm.H
                              Single_Virtual_Correction.H
                              Single_Real_Correction.H)

sherpa_create_git_info(DipoleSubtraction ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(DipoleSubtraction SHARED ${DipoleSubtraction_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(DipoleSubtraction DipoleSubtraction)
target_include_directories(DipoleSubtraction PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
add_dependencies(DipoleSubtraction DipoleSubtraction_Git_Info)
