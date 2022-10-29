
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Cross_Sections/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsXsecs.la)
 set(GITTAG  SHRiMPS/Cross_Sections)
 #include ../../git.make)
set(Cross_Sections_SOURCES  Cross_Sections.C Sigma_Base.C Sigma_Total.C Sigma_Inelastic.C Sigma_Elastic.C Sigma_SD.C Sigma_Partonic.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Cross_Sections)
set(Cross_Sections_HEADERS Cross_Sections.H Sigma_Base.H Sigma_Total.H Sigma_Inelastic.H Sigma_Elastic.H Sigma_SD.H Sigma_Partonic.H)
creategitinfo(Cross_Sections ../Cross_Sections/)
  list(TRANSFORM Cross_Sections_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Cross_Sections/")
 
set(Cross_Sections_esources )
if (SHERPA_BUILD_STATIC)
add_library(ShrimpsXsecs_static STATIC ${ShrimpsXsecs_esources}
                             ${Cross_Sections_SOURCES} 
)
target_include_directories(ShrimpsXsecs_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsXsecs_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsXsecs_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsXsecs_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ShrimpsXsecs)
set_target_properties(ShrimpsXsecs_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(ShrimpsXsecs SHARED ${ShrimpsXsecs_esources}
                             ${Cross_Sections_SOURCES} 
)
target_include_directories(ShrimpsXsecs PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsXsecs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsXsecs DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsXsecs        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsXsecs SOVERSION 0.0.0)
set_target_properties(ShrimpsXsecs        PROPERTIES DEFINE_SYMBOL "")
endif()
