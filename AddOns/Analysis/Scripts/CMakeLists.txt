 
set(Scripts_SOURCES Combine_Analysis.C)
sherpa_create_git_info(Scripts ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_executable(Combine_Analysis  ${Scripts_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
target_include_directories(Combine_Analysis PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Combine_Analysis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_handle_executable(Combine_Analysis)
install(TARGETS Combine_Analysis DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libs)
set_target_properties(Combine_Analysis        PROPERTIES DEFINE_SYMBOL "")
target_link_libraries(Combine_Analysis PRIVATE 
SherpaMain SherpaInitialization SherpaSingleEvents SherpaSoftPhysics SherpaPerturbativePhysics SherpaTools 
ShrimpsMain ShrimpsEvents ShrimpsLadders ShrimpsBeamRemnants ShrimpsXsecs ShrimpsEikonals ShrimpsTools 
HadronsMain HadronsMEs HadronsCurrents HadronsPSs 
Reconnections 
AhadicMain AhadicTools AhadicFormation AhadicDecays 
Amisic AmisicTools AmisicPerturbative 
PhotonsMain PhotonsTools PhotonsPhaseSpace PhotonsMEs PhotonsPhotonSplitter
YFSCEEX YFSMain YFSNLO  YFSTools
Remnants  RemnantsTools 
MCatNLOMain MCatNLOCalculators MCatNLOShowers MCatNLOTools 
DIMTools DIMShower DIMGauge DIMLorentz DIMMain 
CSMain CSCalculators CSShowers CSTools 
DireTools DireShower DireGauge DireLorentz DireMain 
PhasicScales Comix  ComixPhasespace ComixAmplitude 
PhasicMain PhasicChannels PhasicProcess PhasicSelectors PhasicEnhance PhasicDecays
Amegic DipoleSubtraction Amplitude AmegicPSGen 
String Zfunctions 
ExtAmp ExtraXSNLO ExtraXSSpecial ExtraXS2_2 ExtraXS1_2 ExtraXS1_3 ExtraXS 
METoolsExplicit METoolsCurrents METoolsVertices METoolsColors METoolsSpinCorrelations METoolsLoops METoolsMain ModelUFO ModelMain 
PDF 
BeamMain BeamSpectra 
ToolsPhys ToolsMath ToolsYaml ToolsOrg
)
target_link_libraries(Combine_Analysis PRIVATE LibZip::LibZip)
if (SHERPA_ENABLE_HEPMC3)
  target_link_libraries(Combine_Analysis PRIVATE SherpaHepMC3Output  HepMC3::HepMC3 HepMC3::search ${HepMC3_rootIO})
endif()
