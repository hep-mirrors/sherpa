
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#TauPi/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaTauPi.la)
 set(GITTAG  MODEL/TauPi)
 #include ../../git.make)
set(TauPi_SOURCES  TauPi_LC.C Model.C)
creategitinfo(TauPi ../TauPi/)
list(TRANSFORM TauPi_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../TauPi/")
list(APPEND TauPi_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(TauPi_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaTauPi SHARED ${SherpaTauPi_esources}
                             ${TauPi_SOURCES} 
)
target_include_directories(SherpaTauPi PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaTauPi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaTauPi DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaTauPi        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaTauPi SOVERSION 0.0.0)
set_target_properties(SherpaTauPi        PROPERTIES DEFINE_SYMBOL "")
endif()
