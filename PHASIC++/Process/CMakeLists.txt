
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Process/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhasicProcess.la)
 set(GITTAG  PHASIC++/Process)
 #include ../../git.make)
set(Process_SOURCES  Subprocess_Info.C Process_Info.C Process_Base.C Process_Group.C Single_Process.C MCatNLO_Process.C Virtual_ME2_Base.C Tree_ME2_Base.C Color_Correlated_ME2.C Spin_Color_Correlated_ME2.C
 Massive_Kernels.C KP_Terms.C ME_Generators.C ME_Generator_Base.C)
set(Process_HEADERS Subprocess_Info.H Process_Info.H Process_Base.H Process_Group.H Single_Process.H MCatNLO_Process.H Virtual_ME2_Base.H Tree_ME2_Base.H Color_Correlated_ME2.H Spin_Color_Correlated_ME2.H
 External_ME_Args.H Massive_Kernels.H KP_Terms.H ME_Generators.H ME_Generator_Base.H)
creategitinfo(Process ../Process/)
  list(TRANSFORM Process_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Process/")
 
set(Process_esources )
if (SHERPA_BUILD_STATIC)
add_library(PhasicProcess_static STATIC ${PhasicProcess_esources}
                             ${Process_SOURCES} 
)
target_include_directories(PhasicProcess_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicProcess_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicProcess_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicProcess_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PhasicProcess)
set_target_properties(PhasicProcess_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(PhasicProcess SHARED ${PhasicProcess_esources}
                             ${Process_SOURCES} 
)
target_include_directories(PhasicProcess PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhasicProcess PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhasicProcess DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhasicProcess        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhasicProcess SOVERSION 0.0.0)
set_target_properties(PhasicProcess        PROPERTIES DEFINE_SYMBOL "")
endif()
