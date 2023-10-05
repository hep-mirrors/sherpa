set(GRV_SOURCES grvphoton.c
                GRVph_Fortran_Interface.C)

set(GRV_HEADERS GRVph_Fortran_Interface.H)

sherpa_create_git_info(GRV)
 
add_library(GRVSherpa SHARED ${GRV_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_include_directories(GRVSherpa PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../libf2c)
target_link_libraries(GRVSherpa libf2cSherpa)
sherpa_handle_shared_library(GRVSherpa GRVSherpa)
