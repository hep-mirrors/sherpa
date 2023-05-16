set(Colors_SOURCES None_CC.C
                   D_CC.C
                   T_CC.C
                   F_CC.C
                   FF_CC.C
                   TT_CC.C
                   ST_CC.C
                   SF_CC.C
		   SD_CC.C)

set(Colors_HEADERS)

sherpa_create_git_info(Colors)
 
add_library(METoolsColors SHARED ${Colors_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsColors METoolsColors)
