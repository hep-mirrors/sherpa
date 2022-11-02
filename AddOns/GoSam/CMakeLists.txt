
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#GoSam/Makefile.am
 
 set(GITTAG  AddOns/GoSam)
 #include ../../git.make)
set(GoSam_SOURCES  GoSam_Virtual.C GoSam_Interface.C)
set(GOSAM_ADDS   GoSam_Virtual.H GoSam_Interface.H)
 if (SHERPA_ENABLE_GOSAM)
 #pkglib_LTLIBRARIES = libSherpaGoSam.la)
set(GoSam_SOURCES ${GoSam_SOURCES})
 #GOSAM_#EXTRA_DIST = $(GOSAM_ADDS)
 else()
 #GOSAM_#EXTRA_DIST = ${GoSam_SOURCES} $(GOSAM_ADDS)
 endif()
 #EXTRA_DIST = $(#GOSAM_#EXTRA_DIST)
creategitinfo(GoSam ../GoSam/)
list(TRANSFORM GoSam_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../GoSam/")
list(APPEND GoSam_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaGoSam SHARED 
                             ${GoSam_SOURCES} 
)
target_include_directories(SherpaGoSam PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaGoSam DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaGoSam        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaGoSam SOVERSION 0.0.0)
set_target_properties(SherpaGoSam        PROPERTIES DEFINE_SYMBOL "")
