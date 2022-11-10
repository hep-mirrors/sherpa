
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
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
list(APPEND Currents_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsCurrents SHARED 
                             ${Currents_SOURCES} 
)
target_include_directories(METoolsCurrents PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(METoolsCurrents)
install(TARGETS METoolsCurrents DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsCurrents        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsCurrents SOVERSION 0.0.0)
set_target_properties(METoolsCurrents        PROPERTIES DEFINE_SYMBOL "")
