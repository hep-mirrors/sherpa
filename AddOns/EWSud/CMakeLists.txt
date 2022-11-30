##### Automatically or semiautomaticaly generated, do not edit.#########
#EWSud/Makefile.am
 
set(GITTAG  AddOns/EWSud)

set(EWSud_SOURCES Variation_Generator.C
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

set(EWSud_HEADERS Variation_Generator.H
                  Coefficient_Checker.H
                  KFactor_Checker.H
                  Comix_Interface.H
                  HE_Comix_Interface.H
                  EWGroupConstants.H
                  EWSud.H
                  Amplitudes.H
                  KFactor.H
                  Calculator.H
                  Clustered_Calculator.H)

sherpa_create_git_info(EWSud)
 
add_library(SherpaEWSud SHARED ${EWSud_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaEWSud SherpaEWSud)
