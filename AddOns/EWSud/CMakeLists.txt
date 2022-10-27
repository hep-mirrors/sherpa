add_library(SherpaEWSud SHARED 
        Variation_Generator.C 
    Coefficient_Checker.C 
    KFactor_Checker.C 
    Comix_Interface.C 
    HE_Comix_Interface.C 
    EWGroupConstants.C 
    EWSud.C 
    Amplitudes.C  
    KFactor.C 
    Calculator.C 
    Clustered_Calculator.C)

install(TARGETS SherpaEWSud DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaEWSud      PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaEWSud SOVERSION 0.0.0)
