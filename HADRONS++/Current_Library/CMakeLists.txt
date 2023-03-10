set(Current_Library_SOURCES Current_Base.C
                            VA_P_P.C
                            VA_P_S.C
                            VA_P_V.C
                            VA_P_T.C
                            VA_P_A.C
                            VA_0_P.C
                            VA_0_V.C
                            VA_0_PP.C
                            VA_0_PP_strange.C
                            VA_0_PPP.C
                            VA_0_PiPiPiPi3Charged.C
                            VA_0_PiPiPiPi1Charged.C
                            VA_B_DPi.C
                            VA_B_DstarPi.C
                            VA_B_B.C
                            VA_B_B3.C
                            VA_B_B_hepph9409272.C
                            VA_F_F.C)

set(Current_Library_HEADERS Current_Base.H
                            VA_P_P.H
                            VA_P_S.H
                            VA_P_V.H
                            VA_P_T.H
                            VA_P_A.H
                            VA_0_P.H
                            VA_0_V.H
                            VA_0_PP.H
                            VA_0_PP_strange.H
                            VA_0_PPP.H
                            VA_0_PiPiPiPi3Charged.H
                            VA_0_PiPiPiPi1Charged.H
                            VA_B_DPi.H
                            VA_B_DstarPi.H
                            VA_B_B.H
                            VA_B_B3.H
                            VA_F_F.H)

sherpa_create_git_info(Current_Library)
 
add_library(HadronsCurrents SHARED ${Current_Library_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(HadronsCurrents HadronsCurrents)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PhaseSpaceFunctions DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc )
