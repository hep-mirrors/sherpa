
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS.la)
 set(GITTAG  EXTRA_XS/Main)
 #include ../../git.make)
set(Main_SOURCES  Single_Process.C Process_Group.C Simple_XS.C ME_Tools.C ME2_Base.C)
set(Main_HEADERS Single_Process.H Process_Group.H ME_Tools.H ME2_Base.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXS SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(ExtraXS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ExtraXS)
install(TARGETS ExtraXS DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXS SOVERSION 0.0.0)
set_target_properties(ExtraXS        PROPERTIES DEFINE_SYMBOL "")
