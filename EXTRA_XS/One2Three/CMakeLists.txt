
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#One2Three/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS1_3.la)
 set(GITTAG  EXTRA_XS/One2Three)
 #include ../../git.make)
set(One2Three_SOURCES  Comix1to3.C)
set(One2Three_HEADERS Comix1to3.H)
creategitinfo(One2Three ../One2Three/)
list(TRANSFORM One2Three_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../One2Three/")
list(APPEND One2Three_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXS1_3 SHARED 
                             ${One2Three_SOURCES} 
)
target_include_directories(ExtraXS1_3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXS1_3 DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS1_3        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXS1_3 SOVERSION 0.0.0)
set_target_properties(ExtraXS1_3        PROPERTIES DEFINE_SYMBOL "")
