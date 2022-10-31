
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libAmisic.la)
 set(GITTAG  AMISIC/Main)
 #include ../../git.make)
set(Main_SOURCES  Amisic.C)
set(Main_HEADERS Amisic.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Main_esources )
if (SHERPA_BUILD_SHARED)
add_library(Amisic SHARED ${Amisic_esources}
                             ${Main_SOURCES} 
)
target_include_directories(Amisic PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Amisic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS Amisic DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Amisic        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Amisic SOVERSION 0.0.0)
set_target_properties(Amisic        PROPERTIES DEFINE_SYMBOL "")
endif()
