#Manually edited
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#GRS/Makefile.am
 
 #pkglib_LTLIBRARIES = libGRSSherpa.la)
 set(GITTAG  PDF/GRS)
 #include ../../git.make)
 #AM_FFLAGS = -std = legacy)
set(GRS_SOURCES  grsg99.f GRSph_Fortran_Interface.C)
set(GRS_HEADERS GRSph_Fortran_Interface.H)
#nobase_#dist_pkgdata_DATA = GRSGrid/grsg99d.grid GRSGrid/grsg99lo.grid GRSGrid/grsg99m.grid GRSGrid/grsg99nlo.grid
creategitinfo(GRS ../GRS/)
list(TRANSFORM GRS_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../GRS/")
list(APPEND GRS_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(GRS_esources )
if (SHERPA_BUILD_SHARED)
add_library(GRSSherpa SHARED ${GRSSherpa_esources}
                             ${GRS_SOURCES} 
)
target_include_directories(GRSSherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(GRSSherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS GRSSherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(GRSSherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME GRSSherpa SOVERSION 0.0.0)
set_target_properties(GRSSherpa        PROPERTIES DEFINE_SYMBOL "")
target_compile_options(GRSSherpa PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>)
endif()
