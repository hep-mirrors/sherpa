
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#HEFT/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaHEFT.la)
 set(GITTAG  MODEL/HEFT)
 #include ../../git.make)
set(HEFT_SOURCES  Effective_Higgs_Coupling.C Model.C SVV_LC.C SVVV_LC.C SVVP_LC.C SVVVV_LC.C)
set(HEFT_HEADERS Effective_Higgs_Coupling.H)
creategitinfo(HEFT ../HEFT/)
list(TRANSFORM HEFT_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../HEFT/")
list(APPEND HEFT_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaHEFT SHARED 
                             ${HEFT_SOURCES} 
)
target_include_directories(SherpaHEFT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaHEFT DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaHEFT        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaHEFT SOVERSION 0.0.0)
set_target_properties(SherpaHEFT        PROPERTIES DEFINE_SYMBOL "")
