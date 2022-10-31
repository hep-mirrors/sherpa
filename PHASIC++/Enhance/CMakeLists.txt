
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Enhance/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicEnhance.la)
 set(GITTAG  PHASIC++/Enhance)
 #include ../../git.make)
set(Enhance_SOURCES  Enhance_Observable_Base.C Variable_Enhance_Observable.C Fastjet_Enhance_Observable.C)
set(Enhance_HEADERS Enhance_Observable_Base.H)
creategitinfo(Enhance ../Enhance/)
list(TRANSFORM Enhance_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Enhance/")
list(APPEND Enhance_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Enhance_esources )
if (SHERPA_BUILD_SHARED)
add_library(PhasicEnhance SHARED ${PhasicEnhance_esources}
                             ${Enhance_SOURCES} 
)
target_include_directories(PhasicEnhance PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicEnhance PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicEnhance DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicEnhance        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhasicEnhance SOVERSION 0.0.0)
set_target_properties(PhasicEnhance        PROPERTIES DEFINE_SYMBOL "")
endif()
