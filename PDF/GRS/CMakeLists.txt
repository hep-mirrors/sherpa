set(GRS_SOURCES grsg99.c
                GRSph_Fortran_Interface.C)

set(GRS_HEADERS GRSph_Fortran_Interface.H)

sherpa_create_git_info(GRS)
 
add_library(GRSSherpa SHARED ${GRS_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_include_directories(GRSSherpa PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../libf2c)
target_link_libraries(GRSSherpa libf2cSherpa)
sherpa_handle_shared_library(GRSSherpa GRSSherpa)
