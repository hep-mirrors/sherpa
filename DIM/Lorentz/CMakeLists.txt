
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Lorentz/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMLorentz.la)
 set(GITTAG  DIM/Lorentz)
 #include ../../git.make)
set(Lorentz_SOURCES  FFV_FF.C FFV_FI.C FFV_IF.C FFV_II.C VVV_FF.C VVV_FI.C VVV_IF.C VVV_II.C)
creategitinfo(Lorentz ../Lorentz/)
list(TRANSFORM Lorentz_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Lorentz/")
list(APPEND Lorentz_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Lorentz_esources )
if (SHERPA_BUILD_SHARED)
add_library(DIMLorentz SHARED ${DIMLorentz_esources}
                             ${Lorentz_SOURCES} 
)
target_include_directories(DIMLorentz PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DIMLorentz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DIMLorentz DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMLorentz        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMLorentz SOVERSION 0.0.0)
set_target_properties(DIMLorentz        PROPERTIES DEFINE_SYMBOL "")
endif()
