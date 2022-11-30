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

sherpa_create_git_info(NNLO)
 
add_library(NNLOqT SHARED ${NNLO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(NNLOqT NNLOqT)
