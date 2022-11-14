##### Automatically or semiautomaticaly generated, do not edit.#########
#NNLO/Makefile.am
 
 set(GITTAG  AddOns/NNLO)
 #pkglib_LTLIBRARIES = libNNLOqT.la)
set(NNLO_SOURCES  polylog.c Tools.C QT_Selector.C DIS_Selector.C Scale_Setter.C coeffqt.C param.C hard.C DY_KFactor.C higgsfullsm.C H_KFactor.C DISinclu.C DIS_KFactor.C DIS_Scale.C)
#NNLOqT_#EXTRA_DIST = Tools.H QT_Selector.H coeffqt.H param.H hard.H higgsfullsm.H DISinclu.H
 #EXTRA_DIST = $(#NNLOqT_#EXTRA_DIST)
sherpa_create_git_info(NNLO ../NNLO/)
list(TRANSFORM NNLO_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../NNLO/")
list(APPEND NNLO_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(NNLOqT SHARED ${NNLO_SOURCES} )
sherpa_handle_shared_library(NNLOqT NNLOqT)
