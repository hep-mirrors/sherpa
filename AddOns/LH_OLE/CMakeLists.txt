##### Automatically or semiautomaticaly generated, do not edit.#########
#LH_OLE/Makefile.am
 
 set(GITTAG  AddOns/LH_OLE)
set(LH_OLE_SOURCES  LH_OLE_Interface.C LH_OLE_Communicator.C)
 set(LHOLE_ADDS   LH_OLE_Communicator.H)
 if (SHERPA_ENABLE_LHOLE)
 #pkglib_LTLIBRARIES = libSherpaLHOLE.la)
set(LH_OLE_SOURCES ${LH_OLE_SOURCES})
 #LHOLE_#EXTRA_DIST = $(LHOLE_ADDS)
 else()
 #LHOLE_#EXTRA_DIST = ${LH_OLE_SOURCES} $(LHOLE_ADDS)
 endif()
 #EXTRA_DIST = $(#LHOLE_#EXTRA_DIST)
sherpa_create_git_info(LH_OLE ../LH_OLE/)
list(TRANSFORM LH_OLE_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LH_OLE/")
list(APPEND LH_OLE_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaLHOLE SHARED ${LH_OLE_SOURCES} )
sherpa_handle_shared_library(SherpaLHOLE SherpaLHOLE)
