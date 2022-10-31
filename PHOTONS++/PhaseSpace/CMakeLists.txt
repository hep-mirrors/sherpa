
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#PhaseSpace/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhotonsPhaseSpace.la)
 set(GITTAG  PHOTONS++/PhaseSpace)
 #include ../../git.make)
set(PhaseSpace_SOURCES  Avarage_Photon_Number.C Generate_Dipole_Photon_Angle.C Generate_Multipole_Photon_Angle.C)
set(PhaseSpace_HEADERS Avarage_Photon_Number.H Generate_Dipole_Photon_Angle.H Generate_Multipole_Photon_Angle.H)
creategitinfo(PhaseSpace ../PhaseSpace/)
list(TRANSFORM PhaseSpace_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PhaseSpace/")
list(APPEND PhaseSpace_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(PhaseSpace_esources )
if (SHERPA_BUILD_SHARED)
add_library(PhotonsPhaseSpace SHARED ${PhotonsPhaseSpace_esources}
                             ${PhaseSpace_SOURCES} 
)
target_include_directories(PhotonsPhaseSpace PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhotonsPhaseSpace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhotonsPhaseSpace DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhotonsPhaseSpace        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhotonsPhaseSpace SOVERSION 0.0.0)
set_target_properties(PhotonsPhaseSpace        PROPERTIES DEFINE_SYMBOL "")
endif()
