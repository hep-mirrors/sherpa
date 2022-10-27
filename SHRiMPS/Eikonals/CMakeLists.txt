
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
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
 
set(Eikonals_esources )
if (SHERPA_BUILD_STATIC)
add_library(ShrimpsEikonals_static STATIC ${ShrimpsEikonals_esources}
                             ${Eikonals_SOURCES} 
)
target_include_directories(ShrimpsEikonals_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsEikonals_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsEikonals_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsEikonals_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ShrimpsEikonals)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ShrimpsEikonals SHARED ${ShrimpsEikonals_esources}
                             ${Eikonals_SOURCES} 
)
target_include_directories(ShrimpsEikonals PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsEikonals PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsEikonals DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsEikonals        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsEikonals SOVERSION 0.0.0)
endif()
