
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#PS_Library/Makefile.am
 
 #pkglib_LTLIBRARIES = libHadronsPSs.la)
 set(GITTAG  HADRONS++/PS_Library)
 #include ../../git.make)
set(PS_Library_SOURCES  HD_PS_Base.C Three_Body_PSs.C Two_Body_PSs.C Four_Body_PSs.C ResonanceFlavour.C)
set(PS_Library_HEADERS HD_PS_Base.H Three_Body_PSs.H Two_Body_PSs.H Four_Body_PSs.H ResonanceFlavour.H)
creategitinfo(PS_Library ../PS_Library/)
list(TRANSFORM PS_Library_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PS_Library/")
list(APPEND PS_Library_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(HadronsPSs SHARED 
                             ${PS_Library_SOURCES} 
)
target_include_directories(HadronsPSs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(HadronsPSs)
install(TARGETS HadronsPSs DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(HadronsPSs        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME HadronsPSs SOVERSION 0.0.0)
set_target_properties(HadronsPSs        PROPERTIES DEFINE_SYMBOL "")
