
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Showers/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOShowers.la)
 set(GITTAG  MCATNLO/Showers)
 #include ../../git.make)
set(Showers_SOURCES  Splitting_Function_Base.C Splitting_Function_Group.C SF_Coupling.C SF_Lorentz.C Sudakov.C Shower.C Kinematics_Base.C)
set(Showers_HEADERS Splitting_Function_Base.H Splitting_Function_Group.H SF_Coupling.H SF_Lorentz.H SF_Key.H Sudakov.H Shower.H Kinematics_Base.H)
creategitinfo(Showers ../Showers/)
  list(TRANSFORM Showers_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Showers/")
 
set(Showers_esources )
if (SHERPA_BUILD_STATIC)
add_library(MCatNLOShowers_static STATIC ${MCatNLOShowers_esources}
                             ${Showers_SOURCES} 
)
target_include_directories(MCatNLOShowers_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOShowers_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOShowers_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOShowers_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME MCatNLOShowers)
endif()
if (SHERPA_BUILD_SHARED)
add_library(MCatNLOShowers SHARED ${MCatNLOShowers_esources}
                             ${Showers_SOURCES} 
)
target_include_directories(MCatNLOShowers PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOShowers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOShowers DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOShowers        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOShowers SOVERSION 0.0.0)
endif()
