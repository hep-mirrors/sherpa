##### Automatically or semiautomaticaly generated, do not edit.#########
#Higgs/Makefile.am
 
 #AUTOMAKE_OPTIONS = subdir-objects)
 set(GITTAG  AddOns/Higgs)
set(Higgs_SOURCES  dilog.C Wrappers.C Higgs_Tree.C Higgs_Virtual.C Higgs_Interface.C Higgs_Selector.C PS2/CG2.C PS2/C2_0.C PS2/C2_1.C PS2/C2_2.C PS3/CG3.C PS3/C3_0.C PS3/C3_1.C PS3/C3_2.C PS3/C3_3.C PS3/C3_4.C
 PS3/C3_5.C PS3/C3_6.C PS3/C3_7.C PS3/C3_8.C PS3/C3_9.C PS3/C3_10.C PS3/C3_11.C PS3/C3_12.C PS3/C3_13.C PS3/C3_14.C)
#HIGGS_#EXTRA_DIST = gggg.C ggggg.C qqgg_loop.C qqggg_loop.C gggamgam1l.C gggamgam2l.C qqgg_tree.C qqggg_tree.C A_higgs.C A_higgs.h A_spin2.C A_spin2.h A_cont.C A_cont.h Ahiggs.C Ahiggs.h Acont.C Acont.h dilog.h Wrappers.H Higgs_Tree.H Higgs_Virtual.H
 #EXTRA_DIST = $(#HIGGS_#EXTRA_DIST)
sherpa_create_git_info(Higgs ../Higgs/)
list(TRANSFORM Higgs_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Higgs/")
list(APPEND Higgs_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaHiggs SHARED ${Higgs_SOURCES} )
sherpa_handle_shared_library(SherpaHiggs SherpaHiggs)
