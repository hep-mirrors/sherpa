
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#SMDM/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaSMDM.la)
 set(GITTAG  MODEL/SMDM)
 #include ../../git.make)
set(SMDM_SOURCES  Model.C)
creategitinfo(SMDM ../SMDM/)
list(TRANSFORM SMDM_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SMDM/")
list(APPEND SMDM_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(SMDM_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaSMDM SHARED ${SherpaSMDM_esources}
                             ${SMDM_SOURCES} 
)
target_include_directories(SherpaSMDM PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSMDM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSMDM DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSMDM        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSMDM SOVERSION 0.0.0)
set_target_properties(SherpaSMDM        PROPERTIES DEFINE_SYMBOL "")
endif()
