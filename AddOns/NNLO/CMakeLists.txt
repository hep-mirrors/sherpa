set(NNLO_SOURCES polylog.c
                 Tools.C
                 QT_Selector.C
                 DIS_Selector.C
                 Scale_Setter.C
                 coeffqt.C
                 param.C
                 hard.C
                 DY_KFactor.C
                 higgsfullsm.C
                 H_KFactor.C
                 DISinclu.C
                 DIS_KFactor.C
                 DIS_Scale.C)

set(NNLO_HEADERS Tools.H
                 QT_Selector.H
                 coeffqt.H
                 param.H
                 hard.H
                 higgsfullsm.H
                 DISinclu.H)

sherpa_create_git_info(NNLO ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(NNLOqT SHARED ${NNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(NNLOqT NNLOqT)
add_dependencies(NNLOqT NNLO_Git_Info)
