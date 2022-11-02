
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libHadronsMain.la)
 set(GITTAG  HADRONS++/Main)
 #include ../../git.make)
set(Main_SOURCES  Hadron_Decay_Map.C Hadron_Decay_Table.C Hadron_Decay_Channel.C Mixing_Handler.C Tools.C)
set(Main_HEADERS Hadron_Decay_Map.H Hadron_Decay_Table.H Hadron_Decay_Channel.H Mixing_Handler.H Tools.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(HadronsMain SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(HadronsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS HadronsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(HadronsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME HadronsMain SOVERSION 0.0.0)
set_target_properties(HadronsMain        PROPERTIES DEFINE_SYMBOL "")
