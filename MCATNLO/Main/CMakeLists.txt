
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOMain.la)
 set(GITTAG  MCATNLO/Main)
 #include ../../git.make)
set(Main_SOURCES  CS_Cluster_Definitions.C CS_Gamma.C CS_MCatNLO.C)
set(Main_HEADERS CS_Cluster_Definitions.H CS_Gamma.H CS_MCatNLO.H)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(MCatNLOMain_static STATIC ${MCatNLOMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(MCatNLOMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME MCatNLOMain)
set_target_properties(MCatNLOMain_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(MCatNLOMain SHARED ${MCatNLOMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(MCatNLOMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(MCatNLOMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS MCatNLOMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOMain SOVERSION 0.0.0)
set_target_properties(MCatNLOMain        PROPERTIES DEFINE_SYMBOL "")
endif()
