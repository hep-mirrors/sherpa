#Manually edited
#####  Automatically or semiautomaticaly generated, do not edit.#
#Scripts/Makefile.am
 
 if (SHERPA_ENABLE_ANALYSIS)
 #bin_PROGRAMS = Combine_Analysis)
set(Scripts_SOURCES Combine_Analysis.C)
 #Combine_Analysis_LDADD = @ATOOLSLIBS@ @METOOLSLIBS@ @MODELLIBS@ @PDFLIBS@ @CONDITIONAL_LHAPDFLIBS@ @LIBZIP_LDFLAGS@ -lm -ldl)
 endif()
sherpa_create_git_info(Scripts ../Scripts/)
list(TRANSFORM Scripts_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Scripts/")
list(APPEND Scripts_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_executable(Combine_Analysis  ${Scripts_SOURCES} )
target_include_directories(Combine_Analysis PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Combine_Analysis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_handle_executable(Combine_Analysis)
install(TARGETS Combine_Analysis DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libs)
set_target_properties(Combine_Analysis        PROPERTIES DEFINE_SYMBOL "")
target_link_libraries(Combine_Analysis PRIVATE 
 ToolsPhys
 ToolsMath
 ToolsYaml
 ToolsOrg
 METoolsExplicit
 METoolsCurrents
 METoolsVertices
 METoolsColors
 METoolsSpinCorrelations
 METoolsLoops
 METoolsMain
 ModelMain
 ModelUFO
 PDF


)
target_link_libraries(Combine_Analysis PRIVATE LibZip::LibZip)
if (SHERPA_ENABLE_HEPMC3)
  target_link_libraries(Combine_Analysis PRIVATE SherpaHepMC3Output  ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB} )
  if (SHERPA_ENABLE_HEPMC3_ROOT)
    target_link_libraries(Combine_Analysis PRIVATE SherpaHepMC3Output  ${HEPMC3_ROOTIO_LIB} )
  endif()
endif()
if (SHERPA_ENABLE_HEPMC2)
  target_link_libraries(Combine_Analysis PRIVATE SherpaHepMCOutput  ${HEPMC2_LIBRARIES} )
endif()

