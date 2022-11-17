#Manually edited
### Automatically or semiautomaticaly generated, do not edit.
#BlackHat/Makefile.am
 
set(GITTAG  AddOns/BlackHat)
set(BlackHat_SOURCES  BlackHat_Tree.C BlackHat_Virtual.C BlackHat_Interface.C)
set(BLACKHAT_ADDS   BlackHat_Tree.H BlackHat_Virtual.H)
 if (SHERPA_ENABLE_BLACKHAT)
set(BlackHat_SOURCES ${BLACKHAT_SOURCES})
 #BLACKHAT_#EXTRA_DIST = $(BLACKHAT_ADDS)
 else()
 #BLACKHAT_#EXTRA_DIST = ${BLACKHAT_SOURCES} $(BLACKHAT_ADDS)
 endif()
 #EXTRA_DIST = $(#BLACKHAT_#EXTRA_DIST)
sherpa_create_git_info(BlackHat ../BlackHat/)
list(TRANSFORM BlackHat_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../BlackHat/")
list(APPEND BlackHat_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaBlackHat SHARED  ${BlackHat_SOURCES} )
target_link_libraries(SherpaBlackHat PRIVATE ${BLACKHAT_LIBRARY} )
sherpa_handle_shared_library(SherpaBlackHat SherpaBlackHat)
