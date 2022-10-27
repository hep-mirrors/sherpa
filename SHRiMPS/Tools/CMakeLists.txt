
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsTools.la)
 set(GITTAG  SHRiMPS/Tools)
 #include ../../git.make)
set(Tools_SOURCES  DEQ_Solver.C Kernels.C MinBias_Parameters.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Tools)
set(Tools_HEADERS Parameter_Structures.H Special_Functions.H DEQ_Solver.H Kernels.H MinBias_Parameters.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(ShrimpsTools_static STATIC ${ShrimpsTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(ShrimpsTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ShrimpsTools)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ShrimpsTools SHARED ${ShrimpsTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(ShrimpsTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsTools SOVERSION 0.0.0)
endif()
