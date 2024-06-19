set(Main_SOURCES Coulomb.C
                 Define_Dipoles.C
                 FSR.C
                 ISR.C
                 Semi.C
                 YFS_Base.C
                 YFS_Form_Factor.C
                 YFS_Handler.C)

set(Main_HEADERS Coulomb.H
                 Define_Dipoles.H
                 FSR.H
                 ISR.H
                 Semi.H
                 YFS_Base.H
                 YFS_Form_Factor.H
                 YFS_Handler.H)


if (SHERPA_ENABLE_GRIFFIN) 
    list(APPEND Main_SOURCES Griffin_Interface.C)
    list(APPEND Main_HEADERS Griffin_Interface.H)
endif()

sherpa_create_git_info(Main)
 
add_library(YFSMain SHARED ${Main_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSMain YFSMain)
