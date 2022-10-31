
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libDIMMain.la)
 set(GITTAG  DIM/Main)
 #include ../../git.make)
set(Main_SOURCES  Gamma.C MCatNLO.C)
set(Main_HEADERS Gamma.H MCatNLO.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Main_esources )
if (SHERPA_BUILD_SHARED)
add_library(DIMMain SHARED ${DIMMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(DIMMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DIMMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DIMMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DIMMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DIMMain SOVERSION 0.0.0)
set_target_properties(DIMMain        PROPERTIES DEFINE_SYMBOL "")
endif()
