
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libRemnantsTools.la)
 set(GITTAG  REMNANTS/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Kinematics_Generator.C Colour_Generator.C Beam_Decorrelator.C Primordial_KPerp.C Form_Factor.C)
set(Tools_HEADERS Kinematics_Generator.H Colour_Generator.H Beam_Decorrelator.H Primordial_KPerp.H Form_Factor.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(RemnantsTools_static STATIC ${RemnantsTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(RemnantsTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(RemnantsTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS RemnantsTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(RemnantsTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME RemnantsTools)
set_target_properties(RemnantsTools_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(RemnantsTools SHARED ${RemnantsTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(RemnantsTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(RemnantsTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS RemnantsTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(RemnantsTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME RemnantsTools SOVERSION 0.0.0)
set_target_properties(RemnantsTools        PROPERTIES DEFINE_SYMBOL "")
endif()
