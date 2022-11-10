
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND Explicit_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsExplicit SHARED 
                             ${Explicit_SOURCES} 
)
target_include_directories(METoolsExplicit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(METoolsExplicit)
install(TARGETS METoolsExplicit DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsExplicit        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsExplicit SOVERSION 0.0.0)
set_target_properties(METoolsExplicit        PROPERTIES DEFINE_SYMBOL "")
