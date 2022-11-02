#Manually edited
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#OpenLoops/Makefile.am
 
 set(GITTAG  AddOns/OpenLoops)
 #include ../../git.make)
set(OpenLoops_SOURCES  OpenLoops_Virtual.C OpenLoops_Born.C OpenLoops_Interface.C Color_Correlated_ME2.C Spin_Color_Correlated_ME2.C GGH_Process_Manager.C GGH_KFactor_Setter.C GGH_Decay_KFactor_Setter.C)
set(OpenLoops_HEADERS OpenLoops_Virtual.H OpenLoops_Born.H OpenLoops_Interface.H Color_Correlated_ME2.H Spin_Color_Correlated_ME2.H GGH_Process_Manager.H GGH_KFactor_Setter.H)
# libSherpaOpenLoops_la_LDFLAGS = $(AM_LDFLAGS) $(CONDITIONAL_OPENLOOPSLDFLAGS)
 if (SHERPA_SUPPORT_OPENLOOPS)
 #pkglib_LTLIBRARIES = libSherpaOpenLoops.la)
set(OpenLoops_SOURCES ${OpenLoops_SOURCES})
# OPENLOOPS_#EXTRA_DIST = $(OPENLOOPS_ADDS)
 else()
# OPENLOOPS_#EXTRA_DIST = $(OPENLOOPS_SOURCES) $(OPENLOOPS_ADDS)
 endif()
 #EXTRA_DIST = $(OPENLOOPS_#EXTRA_DIST)
creategitinfo(OpenLoops ../OpenLoops/)
list(TRANSFORM OpenLoops_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../OpenLoops/")
list(APPEND OpenLoops_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaOpenLoops SHARED 
                             ${OpenLoops_SOURCES} 
)
target_include_directories(SherpaOpenLoops PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaOpenLoops DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaOpenLoops        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaOpenLoops SOVERSION 0.0.0)
set_target_properties(SherpaOpenLoops        PROPERTIES DEFINE_SYMBOL "")
