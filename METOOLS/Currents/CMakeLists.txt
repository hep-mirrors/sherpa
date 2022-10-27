
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Currents/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsCurrents.la)
 set(GITTAG  METOOLS/Currents)
 #include ../../git.make)
set(Currents_SOURCES  C_Scalar.C C_Spinor.C C_Vector.C C_Pseudo.C S_C.C F_C.C V_C.C P_C.C)
set(Currents_HEADERS C_Scalar.H C_Spinor.H C_Vector.H C_Pseudo.H)
creategitinfo(Currents ../Currents/)
  list(TRANSFORM Currents_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Currents/")
 
set(Currents_esources )
if (SHERPA_BUILD_STATIC)
add_library(METoolsCurrents_static STATIC ${METoolsCurrents_esources}
                             ${Currents_SOURCES} 
)
target_include_directories(METoolsCurrents_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsCurrents_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsCurrents_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsCurrents_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME METoolsCurrents)
endif()
if (SHERPA_BUILD_SHARED)
add_library(METoolsCurrents SHARED ${METoolsCurrents_esources}
                             ${Currents_SOURCES} 
)
target_include_directories(METoolsCurrents PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsCurrents PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsCurrents DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsCurrents        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsCurrents SOVERSION 0.0.0)
endif()
