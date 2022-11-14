##### Automatically or semiautomaticaly generated, do not edit.#########
#GoSam/Makefile.am
 
 set(GITTAG  AddOns/GoSam)
set(GoSam_SOURCES  GoSam_Virtual.C GoSam_Interface.C)
set(GOSAM_ADDS   GoSam_Virtual.H GoSam_Interface.H)
 if (SHERPA_ENABLE_GOSAM)
set(GoSam_SOURCES ${GoSam_SOURCES})
 #GOSAM_#EXTRA_DIST = $(GOSAM_ADDS)
 else()
 #GOSAM_#EXTRA_DIST = ${GoSam_SOURCES} $(GOSAM_ADDS)
 endif()
 #EXTRA_DIST = $(#GOSAM_#EXTRA_DIST)
sherpa_create_git_info(GoSam ../GoSam/)
list(TRANSFORM GoSam_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../GoSam/")
list(APPEND GoSam_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaGoSam SHARED ${GoSam_SOURCES} )
sherpa_handle_shared_library(SherpaGoSam SherpaGoSam)
