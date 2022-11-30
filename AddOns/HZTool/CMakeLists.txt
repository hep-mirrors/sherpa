set(HZTool_SOURCES hznorm.F
                   hztool_sherpa.F
                   HZTool_Interface.C)

set(HZTool_HEADERS HZTool_Wrapper.H)

sherpa_create_git_info(HZTool)
 
add_library(SherpaHZToolAnalysis SHARED ${HZTool_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaHZToolAnalysis SherpaHZToolAnalysis)

target_include_directories( SherpaHZToolAnalysis PRIVATE ${HZTOOL_INCLUDE_DIR} )
target_link_libraries( SherpaHZToolAnalysis PRIVATE ${HZTOOL_LIBRARY} 
${CERN_kernlib_DYNAMIC_LIBRARY}
${CERN_mathlib_DYNAMIC_LIBRARY}
${CERN_packlib_DYNAMIC_LIBRARY}
${MOTIF_LIBRARIES}
)
