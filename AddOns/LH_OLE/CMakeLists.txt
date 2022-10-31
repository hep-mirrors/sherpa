
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#LH_OLE/Makefile.am
 
 set(GITTAG  AddOns/LH_OLE)
 #include ../../git.make)
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
creategitinfo(LH_OLE ../LH_OLE/)
list(TRANSFORM LH_OLE_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LH_OLE/")
list(APPEND LH_OLE_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(LH_OLE_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaLHOLE SHARED ${SherpaLHOLE_esources}
                             ${LH_OLE_SOURCES} 
)
target_include_directories(SherpaLHOLE PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaLHOLE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaLHOLE DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaLHOLE        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaLHOLE SOVERSION 0.0.0)
set_target_properties(SherpaLHOLE        PROPERTIES DEFINE_SYMBOL "")
endif()
