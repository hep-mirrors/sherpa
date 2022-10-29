
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Showers/Makefile.am
 
 #pkglib_LTLIBRARIES = libCSShowers.la)
 set(GITTAG  CSSHOWER++/Showers)
 #include ../../git.make)
set(Showers_SOURCES  Splitting_Function_Base.C Splitting_Function_Group.C SF_Coupling.C SF_Lorentz.C Sudakov.C Shower.C Kinematics_Base.C)
set(Showers_HEADERS Splitting_Function_Base.H Splitting_Function_Group.H SF_Coupling.H SF_Lorentz.H SF_Key.H Sudakov.H Shower.H Kinematics_Base.H)
creategitinfo(Showers ../Showers/)
  list(TRANSFORM Showers_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Showers/")
 
set(Showers_esources )
if (SHERPA_BUILD_STATIC)
add_library(CSShowers_static STATIC ${CSShowers_esources}
                             ${Showers_SOURCES} 
)
target_include_directories(CSShowers_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSShowers_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSShowers_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSShowers_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME CSShowers)
set_target_properties(CSShowers_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(CSShowers SHARED ${CSShowers_esources}
                             ${Showers_SOURCES} 
)
target_include_directories(CSShowers PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CSShowers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CSShowers DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CSShowers        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CSShowers SOVERSION 0.0.0)
set_target_properties(CSShowers        PROPERTIES DEFINE_SYMBOL "")
endif()
