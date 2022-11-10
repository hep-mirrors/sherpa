
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Eikonals/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsEikonals.la)
 set(GITTAG  SHRiMPS/Eikonals)
 #include ../../git.make)
set(Eikonals_SOURCES  Eikonal_Weights.C Rapidity_Density.C Omega_ik.C Single_Channel_Eikonal.C Eikonal_Contributor.C Grid.C Analytic_Eikonal.C Analytic_Contributor.C Eikonal_Creator.C Form_Factors.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Eikonals)
set(Eikonals_HEADERS Eikonal_Weights.H Rapidity_Density.H Omega_ik.H Single_Channel_Eikonal.H Eikonal_Contributor.H Grid.H Analytic_Eikonal.H Analytic_Contributor.H Eikonal_Creator.H Form_Factors.H)
creategitinfo(Eikonals ../Eikonals/)
list(TRANSFORM Eikonals_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Eikonals/")
list(APPEND Eikonals_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsEikonals SHARED 
                             ${Eikonals_SOURCES} 
)
target_include_directories(ShrimpsEikonals PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ShrimpsEikonals)
install(TARGETS ShrimpsEikonals DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsEikonals        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsEikonals SOVERSION 0.0.0)
set_target_properties(ShrimpsEikonals        PROPERTIES DEFINE_SYMBOL "")
