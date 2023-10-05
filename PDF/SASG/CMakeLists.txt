set(SASG_SOURCES sasgam2.c
                 SASGph_Fortran_Interface.C)

set(SASG_HEADERS SASGph_Fortran_Interface.H)

sherpa_create_git_info(SASG)
 
add_library(SASGSherpa SHARED ${SASG_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_include_directories(SASGSherpa PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../libf2c)
target_link_libraries(SASGSherpa libf2cSherpa)
sherpa_handle_shared_library(SASGSherpa SASGSherpa)
