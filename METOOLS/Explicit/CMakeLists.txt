
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Explicit/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsExplicit.la)
 set(GITTAG  METOOLS/Explicit)
 #include ../../git.make)
set(Explicit_SOURCES  C_Object.C Current.C Vertex_Key.C Color_Calculator.C Lorentz_Calculator.C NLO_Counter_Terms.C Dipole_Kinematics.C Dipole_Color.C Dipole_Terms.C Vertex.C)
set(Explicit_HEADERS C_Object.H Current.H Vertex_Key.H Color_Calculator.H Lorentz_Calculator.H NLO_Counter_Terms.H Dipole_Info.H Dipole_Kinematics.H Dipole_Color.H Dipole_Terms.H Vertex.H)
creategitinfo(Explicit ../Explicit/)
  list(TRANSFORM Explicit_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Explicit/")
 
set(Explicit_esources )
if (SHERPA_BUILD_STATIC)
add_library(METoolsExplicit_static STATIC ${METoolsExplicit_esources}
                             ${Explicit_SOURCES} 
)
target_include_directories(METoolsExplicit_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsExplicit_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsExplicit_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsExplicit_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME METoolsExplicit)
endif()
if (SHERPA_BUILD_SHARED)
add_library(METoolsExplicit SHARED ${METoolsExplicit_esources}
                             ${Explicit_SOURCES} 
)
target_include_directories(METoolsExplicit PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(METoolsExplicit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsExplicit DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsExplicit        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsExplicit SOVERSION 0.0.0)
endif()
