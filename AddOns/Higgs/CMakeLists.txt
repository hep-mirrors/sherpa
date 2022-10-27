add_library(SherpaHiggs SHARED 
    dilog.C 
    Wrappers.C 
    Higgs_Tree.C 
    Higgs_Virtual.C 
    Higgs_Interface.C 
    Higgs_Selector.C 
    PS2/CG2.C 
    PS2/C2_0.C 
    PS2/C2_1.C 
    PS2/C2_2.C 
    PS3/CG3.C 
    PS3/C3_0.C 
    PS3/C3_1.C 
    PS3/C3_2.C 
    PS3/C3_3.C 
    PS3/C3_4.C 
    PS3/C3_5.C 
    PS3/C3_6.C 
    PS3/C3_7.C 
    PS3/C3_8.C 
    PS3/C3_9.C 
    PS3/C3_10.C 
    PS3/C3_11.C 
    PS3/C3_12.C 
    PS3/C3_13.C 
    PS3/C3_14.C)

install(TARGETS SherpaHiggs DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaHiggs       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaHiggs SOVERSION 0.0.0)
