set(CT14_SOURCES CT14Pdf.c
                 CT14_Fortran_Interface.C)

sherpa_create_git_info(CT14)
 
add_library(CT14Sherpa SHARED ${CT14_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_include_directories(CT14Sherpa PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../libf2c)
target_link_libraries(CT14Sherpa libf2cSherpa)
sherpa_handle_shared_library(CT14Sherpa CT14Sherpa)
