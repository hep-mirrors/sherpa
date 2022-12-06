set(BlackHat_SOURCES BlackHat_Tree.C
                     BlackHat_Virtual.C
                     BlackHat_Interface.C)

set(BlackHat_HEADERS BlackHat_Tree.H
                     BlackHat_Virtual.H)

sherpa_create_git_info(BlackHat)
 
add_library(SherpaBlackHat SHARED ${BlackHat_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaBlackHat SherpaBlackHat)
target_link_libraries(SherpaBlackHat PRIVATE ${BLACKHAT_LIBRARIES} )
