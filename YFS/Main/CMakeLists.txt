set(YFSMain_SOURCES Coulomb.C
                 Define_Dipoles.C
                 FSR.C
                 ISR.C
                 Semi.C
                 YFS_Base.C
                 YFS_Form_Factor.C
                 YFS_Handler.C)

set(YFSMain_HEADERS Coulomb.H
                 Define_Dipoles.H
                 FSR.H
                 ISR.H
                 Semi.H
                 YFS_Base.H
                 YFS_Form_Factor.H
                 YFS_Handler.H)


sherpa_create_git_info(YFSMain ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(YFSMain SHARED ${YFSMain_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSMain YFSMain)
add_dependencies(YFSMain YFSMain_Git_Info)
# if (SHERPA_ENABLE_LOOPTOOLS) 
#     target_include_directories(YFSMain  PRIVATE ${LOOPTOOLS_DIRS}/include)
#     target_link_libraries(YFSMain  PRIVATE ${LOOPTOOLS_DIRS}/lib/libooptools.a)
# endif()
