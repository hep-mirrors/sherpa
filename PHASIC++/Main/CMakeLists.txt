
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicMain.la)
 set(GITTAG  PHASIC++/Main)
 #include ../../git.make)
set(Main_SOURCES  Color_Integrator.C Helicity_Integrator.C Process_Integrator.C Channel_Creator.C Phase_Space_Point.C Phase_Space_Enhance.C Phase_Space_Handler.C Phase_Space_Integrator.C)
set(Main_HEADERS Color_Integrator.H Helicity_Integrator.H Process_Integrator.H Channel_Creator.H Phase_Space_Point.H Phase_Space_Enhance.H Phase_Space_Variation_Weights.H Phase_Space_Handler.H Phase_Space_Integrator.H)
 #dist_bin_SCRIPTS = plot_stats.sh)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Main_esources )
if (SHERPA_BUILD_SHARED)
add_library(PhasicMain SHARED ${PhasicMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(PhasicMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhasicMain SOVERSION 0.0.0)
set_target_properties(PhasicMain        PROPERTIES DEFINE_SYMBOL "")
endif()
