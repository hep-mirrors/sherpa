
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libComix.la)
 set(GITTAG  COMIX/Main)
 #include ../../git.make)
set(Main_SOURCES  Process_Base.C Process_Group.C Single_Process.C Single_Dipole_Term.C Comix.C)
set(Main_HEADERS Process_Base.H Process_Group.H Single_Process.H Single_Dipole_Term.H)
 #dist_bin_SCRIPTS = init_nlo.sh)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(Comix SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(Comix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(Comix)
install(TARGETS Comix DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Comix        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Comix SOVERSION 0.0.0)
set_target_properties(Comix        PROPERTIES DEFINE_SYMBOL "")
