#Manually edited
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Scripts/Makefile.am
 
 if (SHERPA_ENABLE_ANALYSIS)
 #bin_PROGRAMS = Combine_Analysis)
set(Scripts_SOURCES Combine_Analysis.C)
 #Combine_Analysis_LDADD = @ATOOLSLIBS@ @METOOLSLIBS@ @MODELLIBS@ @PDFLIBS@ @CONDITIONAL_LHAPDFLIBS@ @LIBZIP_LDFLAGS@ -lm -ldl)
 endif()
creategitinfo(Scripts ../Scripts/)
list(TRANSFORM Scripts_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Scripts/")
list(APPEND Scripts_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Scripts_esources )
add_executable(Combine_Analysis  ${Combine_Analysis_esources}
                             ${Scripts_SOURCES} 
)
target_include_directories(Combine_Analysis PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Combine_Analysis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(Combine_Analysis)
install(TARGETS Combine_Analysis DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libs)
set_target_properties(Combine_Analysis        PROPERTIES DEFINE_SYMBOL "")
target_link_libraries(Combine_Analysis PRIVATE 

METoolsColors  METoolsCurrents  METoolsExplicit  METoolsLoops METoolsMain  METoolsSpinCorrelations  METoolsVertices
Amegic Amplitude AmegicPSGen Zfunctions DipoleSubtraction
Remnants RemnantsTools  
PhotonsMain PhotonsMEs PhotonsPhaseSpace PhotonsTools
PhasicChannels PhasicDecays PhasicMain PhasicEnhance PhasicProcess PhasicScales PhasicSelectors
ShrimpsLadders   ShrimpsBeamRemnants  ShrimpsXsecs  ShrimpsEikonals  ShrimpsEvents  ShrimpsMain  ShrimpsTools 
ExtraXS ExtraXS2_2  ExtraXS1_3  ExtraXS1_2
Reconnections  
SherpaMain SherpaTools SherpaSingleEvents SherpaPerturbativePhysics SherpaInitialization LundTools SherpaSoftPhysics
SherpaSM  ModelUFO SherpaTauPi ModelMain
Amisic AmisicTools AmisicPerturbative 
BeamMain BeamSpectra 
Comix ComixAmplitude ComixPhasespace
HadronsMain HadronsCurrents HadronsMEs HadronsPSs
PDF
ToolsMath ToolsOrg ToolsYaml String ToolsPhys
)
target_link_libraries(Combine_Analysis PRIVATE LibZip::LibZip)
if (SHERPA_ENABLE_HEPMC3)
  target_link_libraries(Combine_Analysis PRIVATE SherpaHepMC3Output  ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB} )
  if (SHERPA_ENABLE_HEPMC3)
    target_link_libraries(Combine_Analysis PRIVATE SherpaHepMC3Output  ${HEPMC3_ROOTIO_LIB} )
  endif()
endif()
if (SHERPA_ENABLE_HEPMC2)
  target_link_libraries(Combine_Analysis PRIVATE SherpaHepMCOutput  ${HEPMC2_LIBRARIES} )
endif()

