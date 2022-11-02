
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsMain.la)
 set(GITTAG  SHRiMPS/Main)
 #include ../../git.make)
set(Main_SOURCES  Shrimps.C Cluster_Algorithm.C Soft_Jet_Criterion.C Hadron_Init.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Main)
set(Main_HEADERS Shrimps.H Cluster_Algorithm.H Soft_Jet_Criterion.H Hadron_Init.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsMain SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(ShrimpsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsMain SOVERSION 0.0.0)
set_target_properties(ShrimpsMain        PROPERTIES DEFINE_SYMBOL "")
