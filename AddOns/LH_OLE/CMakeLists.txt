set(LH_OLE_SOURCES LH_OLE_Interface.C
                   LH_OLE_Communicator.C)

set(LH_OLE_HEADERS LH_OLE_Communicator.H)

sherpa_create_git_info(LH_OLE)
 
add_library(SherpaLHOLE SHARED ${LH_OLE_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaLHOLE SherpaLHOLE)
