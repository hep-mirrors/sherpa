set(Main_SOURCES Coulomb.C
                 Define_Dipoles.C
                 FSR.C
                 Griffin_Interface.C
                 ISR.C
                 Semi.C
                 YFS_Base.C
                 YFS_Form_Factor.C
                 YFS_Handler.C)

set(Main_HEADERS Coulomb.H
                 Define_Dipoles.H
                 FSR.H
                 Griffin_Interface.H
                 ISR.H
                 Semi.H
                 YFS_Base.H
                 YFS_Form_Factor.H
                 YFS_Handler.H)


sherpa_create_git_info(Main)
 
add_library(YFSMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSMain YFSMain)
