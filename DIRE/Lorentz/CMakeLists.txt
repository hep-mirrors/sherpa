
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Lorentz/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireLorentz.la)
 set(GITTAG  DIRE/Lorentz)
 #include ../../git.make)
set(Lorentz_SOURCES  FFV_FF.C FFV_FI.C FFV_IF.C FFV_II.C VVV_FF.C VVV_FI.C VVV_IF.C VVV_II.C FFFF_FF.C FFFF_FI.C FFFF_IF.C FFFF_II.C)
creategitinfo(Lorentz ../Lorentz/)
  list(TRANSFORM Lorentz_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Lorentz/")
 
set(Lorentz_esources )
if (SHERPA_BUILD_STATIC)
add_library(DireLorentz_static STATIC ${DireLorentz_esources}
                             ${Lorentz_SOURCES} 
)
target_include_directories(DireLorentz_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireLorentz_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireLorentz_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireLorentz_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME DireLorentz)
set_target_properties(DireLorentz_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(DireLorentz SHARED ${DireLorentz_esources}
                             ${Lorentz_SOURCES} 
)
target_include_directories(DireLorentz PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireLorentz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireLorentz DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireLorentz        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireLorentz SOVERSION 0.0.0)
set_target_properties(DireLorentz        PROPERTIES DEFINE_SYMBOL "")
endif()
