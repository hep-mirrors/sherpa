set(CJK_SOURCES CJKL.c
        cjkhogrid.c
        cjk1grid1.c
        cjk2grid1.c
        CJKph_Fortran_Interface.C)

set(CJK_HEADERS CJKph_Fortran_Interface.H)

sherpa_create_git_info(CJK)

add_library(CJKSherpa SHARED ${CJK_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
target_include_directories(CJKSherpa PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../libf2c)
target_link_libraries(CJKSherpa PRIVATE libf2cSherpa)
sherpa_handle_shared_library(CJKSherpa CJKSherpa)
