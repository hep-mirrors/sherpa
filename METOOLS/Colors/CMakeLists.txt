
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Colors/Makefile.am
 
 #pkglib_LTLIBRARIES = libMEToolsColors.la)
 set(GITTAG  METOOLS/Colors)
 #include ../../git.make)
set(Colors_SOURCES  None_CC.C D_CC.C T_CC.C F_CC.C FF_CC.C TT_CC.C ST_CC.C SF_CC.C)
set(Colors_HEADERS)
creategitinfo(Colors ../Colors/)
list(TRANSFORM Colors_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Colors/")
list(APPEND Colors_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(METoolsColors SHARED 
                             ${Colors_SOURCES} 
)
target_include_directories(METoolsColors PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS METoolsColors DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(METoolsColors        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME METoolsColors SOVERSION 0.0.0)
set_target_properties(METoolsColors        PROPERTIES DEFINE_SYMBOL "")
