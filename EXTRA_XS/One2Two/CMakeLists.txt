
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#One2Two/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS1_2.la)
 set(GITTAG  EXTRA_XS/One2Two)
 #include ../../git.make)
set(One2Two_SOURCES  Comix1to2.C)
set(One2Two_HEADERS Comix1to2.H)
creategitinfo(One2Two ../One2Two/)
list(TRANSFORM One2Two_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../One2Two/")
list(APPEND One2Two_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXS1_2 SHARED 
                             ${One2Two_SOURCES} 
)
target_include_directories(ExtraXS1_2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ExtraXS1_2)
install(TARGETS ExtraXS1_2 DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS1_2        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXS1_2 SOVERSION 0.0.0)
set_target_properties(ExtraXS1_2        PROPERTIES DEFINE_SYMBOL "")
