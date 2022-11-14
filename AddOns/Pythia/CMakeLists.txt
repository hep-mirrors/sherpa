
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Pythia/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaPythia.la)
 if (SHERPA_ENABLE_PYTHIA8)
set(Pythia_SOURCES  Pythia_Jet_Criterion.C Pythia8_Hadronisation.C)
 else()
set(Pythia_SOURCES  Pythia_Jet_Criterion.C)
 #EXTRA_DIST = Pythia8_Hadronisation.C)
 endif()
 set(GITTAG  AddOns/Pythia)
 #include ../../git.make)
creategitinfo(Pythia ../Pythia/)
list(TRANSFORM Pythia_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Pythia/")
list(APPEND Pythia_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaPythia SHARED 
                             ${Pythia_SOURCES} 
)
target_include_directories(SherpaPythia PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SherpaPythia)
install(TARGETS SherpaPythia DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaPythia        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaPythia SOVERSION 0.0.0)
set_target_properties(SherpaPythia        PROPERTIES DEFINE_SYMBOL "")

if (SHERPA_ENABLE_PYTHIA8)   
  target_link_libraries(SherpaPythia PRIVATE ${PYTHIA8_LIBRARIES})
  target_include_directories(SherpaPythia PRIVATE ${PYTHIA8_INCLUDE_DIRS})
endif()
