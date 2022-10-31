#Manually edited
find_package(BlackHat REQUIRED)
message(STATUS "SHERPA: BLACKHAT_LIBRARIES=${BLACKHAT_LIBRARIES} BLACKHAT_INCLUDE_DIRS=${BLACKHAT_INCLUDE_DIRS}")

########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#BlackHat/Makefile.am
 
 set(GITTAG  AddOns/BlackHat)
 #include ../../git.make)
set(BlackHat_SOURCES  BlackHat_Tree.C BlackHat_Virtual.C BlackHat_Interface.C)
set(BLACKHAT_ADDS   BlackHat_Tree.H BlackHat_Virtual.H)
 if (BLACKHAT_SUPPORT)
 #pkglib_LTLIBRARIES = libSherpaBlackHat.la)
set(BlackHat_SOURCES ${BLACKHAT_SOURCES})
 #BLACKHAT_#EXTRA_DIST = $(BLACKHAT_ADDS)
 else()
 #BLACKHAT_#EXTRA_DIST = ${BLACKHAT_SOURCES} $(BLACKHAT_ADDS)
 endif()
 #EXTRA_DIST = $(#BLACKHAT_#EXTRA_DIST)
creategitinfo(BlackHat ../BlackHat/)
list(TRANSFORM BlackHat_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../BlackHat/")
list(APPEND BlackHat_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(BlackHat_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaBlackHat SHARED ${SherpaBlackHat_esources}
                             ${BlackHat_SOURCES} 
)
target_include_directories(SherpaBlackHat PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaBlackHat PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaBlackHat DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaBlackHat        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaBlackHat SOVERSION 0.0.0)
set_target_properties(SherpaBlackHat        PROPERTIES DEFINE_SYMBOL "")
target_link_libraries(SherpaBlackHat PRIVATE ${BLACKHAT_LIBRARIES})
target_include_directories(SherpaBlackHat PRIVATE ${BLACKHAT_INCLUDE_DIRS})
endif()
