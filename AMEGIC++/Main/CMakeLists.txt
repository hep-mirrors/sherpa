
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libAmegic.la)
 set(GITTAG  AMEGIC++/Main)
 #include ../../git.make)
set(Main_SOURCES  ColorSC.C Helicity.C Topology.C Polarisation.C Pol_Info.C Point.C Process_Tags.C Amegic_Base.C Process_Base.C Single_Process.C Single_Process_MHV.C Single_Process_External.C Single_Process_Combined.C
 Process_Group.C Amegic.C)
set(Main_HEADERS Tools.H ColorSC.H Helicity.H Topology.H Polarisation.H Pol_Info.H Point.H Process_Tags.H Amegic_Base.H Process_Base.H Single_Process.H Single_Process_MHV.H Single_Process_External.H Single_Process_Combined.H
 Process_Group.H)
 #dist_pkgdata_SCRIPTS = makelibs)
 #dist_bin_SCRIPTS = make2scons)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Main_esources )
if (SHERPA_BUILD_SHARED)
add_library(Amegic SHARED ${Amegic_esources}
                             ${Main_SOURCES} 
)
target_include_directories(Amegic PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Amegic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS Amegic DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Amegic        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Amegic SOVERSION 0.0.0)
set_target_properties(Amegic        PROPERTIES DEFINE_SYMBOL "")
endif()
