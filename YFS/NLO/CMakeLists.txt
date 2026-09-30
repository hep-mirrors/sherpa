set(YFSNLO_SOURCES NLO_Base.C
                 NLO_Real.C
                 NLO_Virtual.C
                 NLO_RealVirtual.C
                 NLO_RealReal.C
                 NLO_Kinematics.C
                 NLO_LoopHelicity.C
                 NLO_Checks.C
                 Photon_Counterterm.C
                 ME_Compare.C
                 Real_Correction.C
                 RealVirtual.C
                 Virtual.C
                 VirtualVirtual.C)

set(YFSNLO_HEADERS  NLO_Base.H
                 ME_Compare.H
                 Photon_Counterterm.H
                 Real_Correction.H
                 RealVirtual.H
                 Virtual.H
                 VirtualVirtual.H)


sherpa_create_git_info(YFSNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(YFSNLO ${SHERPA_SUBLIB_TYPE} ${YFSNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSNLO YFSNLO)
add_dependencies(YFSNLO YFSNLO_Git_Info)
