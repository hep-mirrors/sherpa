
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#EXTAMP/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtAmp.la)
 set(GITTAG  EXTAMP)
 #include ../git.make)
set(EXTAMP_SOURCES  External_ME_Interface.C Process.C Born_Process.C RS_Process.C BVI_Process.C CS_Dipole.C CS_Dipoles.C Dipole_Wrapper_Process.C)
set(EXTAMP_HEADERS External_ME_Interface.H Process_Group.H Process.H Born_Process.H RS_Process.H BVI_Process.H CS_Dipole.H CS_Dipoles.H Dipole_Wrapper_Process.H)
creategitinfo(EXTAMP ../EXTAMP/)
  list(TRANSFORM EXTAMP_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../EXTAMP/")
 
set(EXTAMP_esources )
if (SHERPA_BUILD_STATIC)
add_library(ExtAmp_static STATIC ${ExtAmp_esources}
                             ${EXTAMP_SOURCES} 
)
target_include_directories(ExtAmp_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtAmp_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtAmp_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtAmp_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ExtAmp)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ExtAmp SHARED ${ExtAmp_esources}
                             ${EXTAMP_SOURCES} 
)
target_include_directories(ExtAmp PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtAmp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtAmp DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtAmp        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtAmp SOVERSION 0.0.0)
endif()
