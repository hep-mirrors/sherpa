cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)
###########################Package@#####################################
project(SHERPA-MC LANGUAGES C  CXX Fortran )
#The complicated way to set the version is used to assure the simple 
#replacement of version variable with sed will not break the CMakeLists.txt
if (NOT __CL_VERSION__)
  set(EXTERNAL_VERSION_STRING "3.0.0alpha1")
else()
  set(EXTERNAL_VERSION_STRING "__CL_VERSION__.1")
endif()
string(REPLACE "."  ";" EXTERNAL_VERSION_LIST "${EXTERNAL_VERSION_STRING}")
list(GET EXTERNAL_VERSION_LIST 0 SHERPA_VERSION_MAJOR)
list(GET EXTERNAL_VERSION_LIST 1 SHERPA_VERSION_MINOR)
list(GET EXTERNAL_VERSION_LIST 2 SHERPA_VERSION_PATCH)
set(PROJECT_VERSION ${SHERPA_VERSION_MAJOR}.${SHERPA_VERSION_MINOR}.${SHERPA_VERSION_PATCH})
SET(SHERPA_VERSION ${SHERPA_VERSION_MAJOR}.${SHERPA_VERSION_MINOR}.${SHERPA_VERSION_PATCH})
message(STATUS "SHERPA: SHERPA_VERSION=${SHERPA_VERSION}")
message(STATUS "SHERPA: SHERPA_VERSION_MAJOR=${SHERPA_VERSION_MAJOR}")
message(STATUS "SHERPA: SHERPA_VERSION_MINOR=${SHERPA_VERSION_MINOR}")
message(STATUS "SHERPA: SHERPA_VERSION_PATCH=${SHERPA_VERSION_PATCH}")
message(STATUS "SHERPA: Building on CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}")
if(CMAKE_BUILD_TYPE)
  message(STATUS "SHERPA: Build type: ${CMAKE_BUILD_TYPE}")
else(CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
  message(STATUS "SHERPA: Build type not specified, use: ${CMAKE_BUILD_TYPE}")
endif(CMAKE_BUILD_TYPE)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

if (APPLE)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
###########################Options#####################################
option(SHERPA_ENABLE_MPI                "Enables building of tests." OFF)
option(SHERPA_ENABLE_HZTOOL                "Enables building of tests." ON)
option(SHERPA_ENABLE_RECOLA                "Enables building of tests." ON)
option(SHERPA_ENABLE_OPENLOOPS                "Enables building of tests." ON)
option(SHERPA_ENABLE_LHOLE                "Enables building of tests." ON)
option(SHERPA_ENABLE_MADLOOP                "Enables building of tests." OFF)
option(SHERPA_ENABLE_BLACKHAT                "Enables building of tests." ON)
option(SHERPA_ENABLE_ROOT                "Enables building of tests." ON)
option(SHERPA_ENABLE_HEPMC3                "Enables building of tests." ON)
option(SHERPA_ENABLE_HEPMC2                "Enables building of tests." ON)
option(SHERPA_ENABLE_PYTHIA                "Enables building of tests." ON)
option(SHERPA_ENABLE_RIVET                "Enables building of tests." ON)
option(SHERPA_ENABLE_LHAPDF                "Enables building of tests." ON)
option(SHERPA_ENABLE_PYTHON                "Enables building of tests." ON)
option(SHERPA_POSITION_INDEPENDENT_CODE  "Position independent code for static libraries." OFF)
option(SHERPA_BUILD_SHARED               "Build shared libraries." ON)
option(SHERPA_BUILD_STATIC               "Build static libraries." OFF)
option(PYTHIA_SUPPORT               "Pythia support" ON)
option(GZIP_SUPPORT                "GZIP support" ON)
###########################Options handeling#####################################
###########Respect HFS, e.g. defines lib or lib64 when it is needed.####
include("GNUInstallDirs")
if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  if (NOT CMAKE_MACOSX_RPATH) 
    set(CMAKE_MACOSX_RPATH 1)
  endif()
  message(STATUS "SHERPA: CMAKE_MACOSX_RPATH  ${CMAKE_MACOSX_RPATH}. Set by default to 1. Set to 0 to disable RPATH in the build executables.")
endif()
#### Modules ###########################################################
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
#### RPATH handling ####################################################
set(CMAKE_MACOSX_RPATH 1)
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")
###################Windows##############################################
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
###################FLAGS################################################
SET(CMAKE_POSITION_INDEPENDENT_CODE ${OLD_CMAKE_POSITION_INDEPENDENT_CODE})
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message(STATUS "SHERPA: Fortran_COMPILER_NAME=${Fortran_COMPILER_NAME}  CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}  CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
message(STATUS "SHERPA: C_COMPILER_NAME=${C_COMPILER_NAME}  CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}   CMAKE_C_COMPILER_VERSION=${CMAKE_C_COMPILER_VERSION}")

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
endif()

message(STATUS "SHERPA: CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "SHERPA: CMAKE_MODULE_LINKER_FLAGS=${CMAKE_MODULE_LINKER_FLAGS}")
message(STATUS "SHERPA: CMAKE_SHARED_LINKER_FLAGS=${CMAKE_SHARED_LINKER_FLAGS}")
message(STATUS "SHERPA: CMAKE_STATIC_LINKER_FLAGS=${CMAKE_STATIC_LINKER_FLAGS}")
message(STATUS "SHERPA: CMAKE_EXE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_EXE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}")
message(STATUS "SHERPA: CMAKE_MODULE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_MODULE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}")
message(STATUS "SHERPA: CMAKE_SHARED_LINKER_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_SHARED_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}")
message(STATUS "SHERPA: CMAKE_STATIC_LINKER_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_STATIC_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}")
message(STATUS "SHERPA: CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
message(STATUS "SHERPA: CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
message(STATUS "SHERPA: CPP=${CPP}")

include(CheckFortranCompilerFlag)
include(CheckCCompilerFlag)
set ( NEW_CMAKE_C_FLAGS " ")
string(REPLACE " " ";" LIST_CMAKE_C_FLAGS "${C_C_FLAGS_TO_CHECK} ${CMAKE_C_FLAGS}")
foreach(fx IN LISTS LIST_CMAKE_C_FLAGS)
  if (fx)
    check_c_compiler_flag(${fx} C_COMPILER_SUPPORTS_${fx})
    if(C_COMPILER_SUPPORTS_${fx})
      list( APPEND NEW_CMAKE_C_FLAGS ${fx})
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${fx}" ) #neede for dependent flags
    endif()
  endif()
endforeach()
list(JOIN NEW_CMAKE_C_FLAGS " " CMAKE_C_FLAGS_X )
if (NOT CMAKE_FORCE_FLAGS)
  set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS_X})
endif()
message(STATUS "SHERPA: CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
########################################################################
set ( NEW_CMAKE_Fortran_FLAGS " ")
string(REPLACE " " ";" LIST_CMAKE_Fortran_FLAGS  "${F_Fortran_FLAGS_TO_CHECK} ${CMAKE_Fortran_FLAGS}")
foreach(fl IN LISTS LIST_CMAKE_Fortran_FLAGS)
  if (fl)
    check_Fortran_compiler_flag(${fl} Fortran_COMPILER_SUPPORTS_${fl})
    if(Fortran_COMPILER_SUPPORTS_${fl})
      list( APPEND NEW_CMAKE_Fortran_FLAGS ${fl})
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${fl}" ) #neede for dependent flags      
    endif()
  endif()
endforeach()
list(JOIN NEW_CMAKE_Fortran_FLAGS " " CMAKE_Fortran_FLAGS_X )
if (NOT CMAKE_FORCE_FLAGS)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_X} -std=legacy")
endif()
message(STATUS "SHERPA : CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
########################################################################
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if(SHERPA_ENABLE_MPI)
 find_package(MPI REQUIRED)
 set (LOADMPIFORPY "from mpi4py import MPI")
else()
 set (LOADMPIFORPY "")
endif()

add_subdirectory(AHADIC++)
add_subdirectory(AMEGIC++)
add_subdirectory(AMISIC++)

add_subdirectory(BEAM)
add_subdirectory(COMIX)
add_subdirectory(CSSHOWER++)
add_subdirectory(DIM)
add_subdirectory(DIRE)
add_subdirectory(EXTAMP)
add_subdirectory(EXTRA_XS)
add_subdirectory(HADRONS++)
add_subdirectory(METOOLS)
add_subdirectory(MCATNLO)
add_subdirectory(MODEL)
add_subdirectory(PHASIC++)
add_subdirectory(PHOTONS++)
add_subdirectory(PDF/CT14  )
add_subdirectory(PDF/Electron)  
add_subdirectory(PDF/GRS  )
add_subdirectory(PDF/GRV  )
IF (SHERPA_ENABLE_LHAPDF)
  find_package(LHAPDF REQUIRED)
  message(STATUS "SHERPA-MC: LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS}")
  add_subdirectory(PDF/LHAPDF)
endif()
add_subdirectory(ATOOLS)
add_subdirectory(PDF/Main)
add_subdirectory(PDF/NNPDF)  
add_subdirectory(PDF/SAL) 
add_subdirectory(PDF/SASG)
add_subdirectory(RECONNECTIONS)
add_subdirectory(REMNANTS)
add_subdirectory(SHRiMPS)
add_subdirectory(SHERPA)

if (SHERPA_ENABLE_HEPMC3 OR SHERPA_ENABLE_HEPMC2)
  if (SHERPA_ENABLE_HEPMC3)
     find_package(HepMC2 REQUIRED)
     message(STATUS "SHERPA: HEPMC2_VERSION=${HEPMC2_VERSION} HEPMC2_LIBRARIES=${HEPMC2_LIBRARIES} HEPMC2_FIO_LIBRARIES=${HEPMC2_FIO_LIBRARIES} HEPMC2_INCLUDE_DIRS=${HEPMC2_INCLUDE_DIRS}")
   endif()
   if (SHERPA_ENABLE_HEPMC3)
    find_package(HepMC3 REQUIRED)
    message(STATUS "SHERPA: HEPMC3_VERSION=${HEPMC3_VERSION} HEPMC3_LIBRARIES=${HEPMC3_LIBRARIES}  HEPMC3_INCLUDE_DIR=${HEPMC3_INCLUDE_DIR}") 
  endif()
  add_subdirectory(AddOns/HepMC)
endif()


set( IDIRS  AHADIC++
AMEGIC++
AMISIC++
ATOOLS
BEAM
COMIX
CSSHOWER++
DIM
DIRE
EXTAMP
EXTRA_XS
HADRONS++
METOOLS
MCATNLO
MODEL
PHASIC++
PHOTONS++
PDF
RECONNECTIONS
REMNANTS
SHRiMPS
SHERPA
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/Python 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC/AddOns/ COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" 
                             PATTERN "*.deps*" EXCLUDE   
                             PATTERN "*.lib*" EXCLUDE 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE 
                             )

  foreach ( IDIR ${IDIRS} )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${IDIR} 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC  COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" PATTERN "*.icc"  
                             PATTERN "*/fjcore.hh" 
                             PATTERN "*/Node.C" 
                             PATTERN "*/Getter_Function.C" 
                             PATTERN "*/Dummy.C" 
                             PATTERN "*/Math/Cluster_Algorithm.C"
                             PATTERN "*.deps*" EXCLUDE   
                             PATTERN "*.lib*" EXCLUDE 
                             PATTERN "*/SoftPhysics/Cluster_Algorithm.H" EXCLUDE
                             PATTERN "*Cluster_Decay_Analysis.H" EXCLUDE 
                             PATTERN "*Splitting_Tools.H" EXCLUDE 
                             PATTERN "*Running_Alpha_HV.H" EXCLUDE 
                             PATTERN "*Sigma_DD.H" EXCLUDE 
                             PATTERN "*Full_Amplitude_Base.H" EXCLUDE 
                             PATTERN "*Event_Generation/Ladder.H" EXCLUDE 
                             PATTERN "*Ladder_Generator_LDC.H" EXCLUDE 
                             PATTERN "*Quasi_Elastic_Event_Generator.H" EXCLUDE 
                             PATTERN "*Double_Diffractive_Event_Generator.H" EXCLUDE 
                             PATTERN "*EXTRA_XS/Special" EXCLUDE 
                             PATTERN "*EXTRA_XS/Two2Two" EXCLUDE 
                             PATTERN "*Math/Sobol/*" EXCLUDE 
                             PATTERN "*Math/Sobol" EXCLUDE 
                             PATTERN "*EXTAMP/Example*" EXCLUDE 
                             PATTERN "*HADRONS++/Run/Main*" EXCLUDE 
                             PATTERN "*MODEL/AQGC*" EXCLUDE 
                             PATTERN "*MODEL/MUED*" EXCLUDE 
                             PATTERN "*MODEL/SM/*" EXCLUDE 
                             PATTERN "*MODEL/SM" EXCLUDE 
                             PATTERN "*MODEL/SMDM*" EXCLUDE 
                             PATTERN "*MODEL/SMGold*" EXCLUDE 
                             PATTERN "*MODEL/TauPi*" EXCLUDE 
                             PATTERN "*DIM/Gauge*" EXCLUDE 
                             PATTERN "*DIM/Lorentz*" EXCLUDE 
                             PATTERN "*DIRE/Gauge*" EXCLUDE 
                             PATTERN "*DIRE/Lorentz*" EXCLUDE 
                             PATTERN "*DIRE/Main*" EXCLUDE 
                             PATTERN "*PDF/CT14*" EXCLUDE 
                             PATTERN "*PDF/GRS/GRSGrid*" EXCLUDE 
                             PATTERN "*PDF/SAL/SALGrid*" EXCLUDE 
                             PATTERN "*PDF/LHAPDF" EXCLUDE 
                             PATTERN "*METOOLS/Vertices" EXCLUDE 
                             PATTERN "*MCATNLO/Calculators" EXCLUDE 
                             PATTERN "*CSSHOWER++/Calculators" EXCLUDE 
                             PATTERN "*METOOLS/Colors" EXCLUDE 
                             PATTERN "*SHERPA/Run" EXCLUDE 
                             PATTERN "*HADRONS++/Run" EXCLUDE 
                             PATTERN "*makefile*" 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE )    
endforeach (  IDIR ${IDIRS}  )

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Examples DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC 
  COMPONENT doc 
                         PATTERN test.py.in EXCLUDE 
                         PATTERN CI EXCLUDE 
                         PATTERN Userhook EXCLUDE 
                         PATTERN HPC EXCLUDE 
                         PATTERN NLO_EW_* EXCLUDE 
                         PATTERN NLO_QCD_* EXCLUDE 
                         PATTERN RemnantTests EXCLUDE 
                         PATTERN VBF_Interference EXCLUDE 
                         PATTERN */FixedOrder_NNLO EXCLUDE 
                         PATTERN */Instanton/LHC_13 EXCLUDE 
                         PATTERN */LHC_HInt_MCatNLO EXCLUDE 
                         PATTERN */Instanton EXCLUDE 
                         PATTERN */PETRA_Jets EXCLUDE 
                         PATTERN */RelicDensity EXCLUDE 
                         PATTERN */SM_ggV EXCLUDE 
                         PATTERN */*UN2LOPS EXCLUDE 
                         PATTERN Examples/Makefile EXCLUDE 
                         PATTERN Makefile.am EXCLUDE 
                         PATTERN Makefile.in EXCLUDE 
                         PATTERN Makefile.in EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.197.yaml  EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.206.yaml  EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.EPA_PDF.yaml  EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.EPA.yaml EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.GammaPDF.yaml EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.PDF.yaml EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.Reconnect.yaml EXCLUDE 
                         PATTERN */Jets_at_LeptonColliders/LEP_Jets/Sherpa.test.yaml EXCLUDE 
                         PATTERN */SingleTop_Channels/Sherpa.tj-s_channel-nf4.yaml EXCLUDE 
                         PATTERN */Soft_QCD/LHC_7TeV_MinBias/Elastic.yaml EXCLUDE 
                         PATTERN */Soft_QCD/LHC_7TeV_MinBias/Shrimps.yaml EXCLUDE
                         PATTERN */Soft_QCD/LHC_7TeV_MinBias/SoftDiffractive.yaml EXCLUDE 
                         PATTERN */Tops_plus_Jets/LHC_Tops/Sherpa.Reconnect.yaml EXCLUDE 
                         PATTERN */VV_plus_Jets/LHC_2l2nuJets_GluonInitiated/Sherpa.tree.yaml EXCLUDE 
                         PATTERN */VV_plus_Jets/LHC_4lJets_GluonInitiated/Sherpa.tree.yaml EXCLUDE 
                          )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CT14/CT14Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc PATTERN "*ct14ll*" EXCLUDE  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/GRS/GRSGrid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/SAL/SALGrid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/HEPEVT DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ATOOLS/Math/Sobol/21201.gz DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC/Sobol COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/COMIX/Main/init_nlo.sh DESTINATION  ${CMAKE_INSTALL_BINDIR} )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/get_bibtex DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/plot_graphs DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2 DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
add_custom_target( unZip ALL)
add_custom_command(TARGET unZip PRE_BUILD
       COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/NNPDF31_nnlo_as_0118_mc
       COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2
       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2
       COMMENT "UnpackingNNPDF31_nnlo_as_0118_mc.tar.bz2"
       VERBATIM)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/NNPDF31_nnlo_as_0118_mc DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc )

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/HADRONS++/Decaydata.zip DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/makelibs   ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/MadLoop/sconsloops DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PHASIC++/Main/plot_stats.sh  
                ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/make2scons   
                ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/Run/Sherpa-config
                DESTINATION  ${CMAKE_INSTALL_BINDIR} )

set(Sherpa_SRCS_C Sherpa.C Filter.C)
list(TRANSFORM Sherpa_SRCS_C PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Main/")
add_executable(Sherpa  ${Sherpa_SRCS_C} ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA//Run/Main.C)
target_link_libraries(Sherpa PRIVATE    
METoolsColors  METoolsCurrents  METoolsExplicit  METoolsLoops METoolsMain  METoolsSpinCorrelations  METoolsVertices
Amegic Amplitude AmegicPSGen Zfunctions DipoleSubtraction
Remnants RemnantsTools  
PhotonsMain PhotonsMEs PhotonsPhaseSpace PhotonsTools
PhasicChannels PhasicDecays PhasicMain PhasicEnhance PhasicProcess PhasicScales PhasicSelectors
ShrimpsLadders   ShrimpsBeamRemnants  ShrimpsXsecs  ShrimpsEikonals  ShrimpsEvents  ShrimpsMain  ShrimpsTools 
ExtraXS ExtraXS2_2  ExtraXS1_3  ExtraXS1_2
Reconnections  
zip   
SherpaMain SherpaTools SherpaSingleEvents SherpaPerturbativePhysics SherpaInitialization LundTools SherpaSoftPhysics
SherpaSM  ModelUFO SherpaTauPi ModelMain
Amisic AmisicTools AmisicPerturbative 
BeamMain BeamSpectra 
Comix ComixAmplitude ComixPhasespace
HadronsMain HadronsCurrents HadronsMEs HadronsPSs
PDF
ToolsMath ToolsOrg ToolsYaml String ToolsPhys
)
if (SHERPA_ENABLE_HEPMC3)
  target_link_libraries(Sherpa PRIVATE SherpaHepMC3Output  ${HEPMC3_LIBRARIES} )
endif()
if (SHERPA_ENABLE_HEPMC2)
  target_link_libraries(Sherpa PRIVATE SherpaHepMCOutput  ${HEPMC2_LIBRARIES} )
endif()

target_include_directories(Sherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
install(TARGETS Sherpa DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libs)

add_subdirectory(AddOns/Weights)
add_subdirectory(AddOns/NNLO)
add_subdirectory(AddOns/Higgs)
if (SHERPA_ENABLE_PYTHON)
  if (NOT  SHERPA-MC_PYTHON_VERSION)
    find_package (Python 3 COMPONENTS  Development Interpreter)
  else()
    find_package (Python ${SHERPA-MC_PYTHON_VERSION} EXACT COMPONENTS  Development Interpreter)
  endif()
  if (NOT SHERPA-MC_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
    set (SHERPA-MC_Python_SITEARCH ${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/)
    message(STATUS "SHERPA-MC python: SHERPA-MC_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR} defaults to ${SHERPA-MC_Python_SITEARCH}")
  else()
    message(STATUS "SHERPA-MC python: SHERPA-MC_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR} is set to to ${SHERPA-MC_Python_SITEARCH}")
  endif()
  add_subdirectory(AddOns/Python)
endif()

if (SHERPA_ENABLE_RIVET)
  find_package(Rivet REQUIRED)
  message(STATUS "SHERPA-MC: RIVET_FOUND=${RIVET_FOUND} RIVET_VERSION=${RIVET_VERSION} RIVET_LIBRARIES=${RIVET_LIBRARIES} RIVET_INCLUDE_DIRS=${RIVET_INCLUDE_DIRS} RIVET_DATA_PATH=${RIVET_DATA_PATH} RIVET_ANALYSIS_PATH=${RIVET_ANALYSIS_PATH} Rivet_HEPMC2_FOUND=${Rivet_HEPMC2_FOUND}  Rivet_HEPMC3_FOUND=${Rivet_HEPMC3_FOUND}")
  add_subdirectory(AddOns/Rivet)
endif()
if (SHERPA_ENABLE_PYTHIA)
  add_subdirectory(AddOns/Pythia)
endif()
if (SHERPA_ENABLE_ROOT)
  add_subdirectory(AddOns/Root)
endif()
if (SHERPA_ENABLE_BLACKHAT)
  add_subdirectory(AddOns/BlackHat)
endif()
if (SHERPA_ENABLE_MADLOOP)
  add_subdirectory(AddOns/MadLoop)
endif()
if (SHERPA_ENABLE_LHOLE)
  add_subdirectory(AddOns/LH_OLE)
endif()
if (SHERPA_ENABLE_OPENLOOPS)
  add_subdirectory(AddOns/OpenLoops)
endif()
if (SHERPA_ENABLE_RECOLA)
  add_subdirectory(AddOns/Recola)
endif()
if (SHERPA_ENABLE_HZTOOL)
  add_subdirectory(AddOns/HZTool)
endif()

########################################################################
# Config files
set (CXX ${CMAKE_CXX_COMPILER})
set (CXXFLAGS ${CMAKE_CXX_FLAGS})
set (CC ${CMAKE_C_COMPILER})
set (CFLAGS ${CMAKE_C_FLAGS})
set (FC ${CMAKE_Fortran_COMPILER})
set (INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC)
set (LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC)
set (DATADIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR})
set (prefix ${CMAKE_INSTALL_PREFIX})
set (PACKAGE_VERSION ${SHERPA_VERSION})
set (host "somehost")

if (NOT HEPEVT_CB_SIZE)
  set(HEPEVT_CB_SIZE  4000)
endif()
set( LD_LIBRARY_PATH "LD_LIBRARY_PATH")
if (Darwin)
  set( LD_LIBRARY_PATH "DYLD_LIBRARY_PATH")
endif()
set (PYLIBDIR ${SHERPA-MC_Python_SITEARCH})
string(REPLACE ";"  " " CONDITIONAL_HEPMC2LIBS " ${HEPMC2_LIBRARIES} ")
string(REPLACE ";"  " " CONDITIONAL_HEPMC3LIBS " ${HEPMC3_LIBRARIES} ")
string(REPLACE ";"  " " CONDITIONAL_LHAPDFLIBS  " ${LHAPDF_LIBRARIES} ")
string(REPLACE ";"  " " CONDITIONAL_RIVETLDADD " ${RIVET_LIBRARIES} ")
string(REPLACE ";"  "  -I" CONDITIONAL_HEPMC2INCS " -I${HEPMC2_INCLUDE_DIRS} ")
string(REPLACE ";"  " -I " CONDITIONAL_HEPMC3INCS " -I${HEPMC3_INCLUDE_DIRS} ")
string(REPLACE ";"  " -I" CONDITIONAL_LHAPDFINCS  " -I${LHAPDF_INCLUDE_DIRS} ")
string(REPLACE ";"  " " CONDITIONAL_RIVETCPPFLAGS "  ")
#  --version)
#    echo @PACKAGE_VERSION@;;
#  --prefix)
#    echo @prefix@;;
#  --libs)
#    echo -L${SHERPA_LIBRARY_PATH:-@LIBDIR@} -lSherpaMain 
#@CONDITIONAL_HEPMC2LIBS@ @CONDITIONAL_HEPMC3LIBS@ @CONDITIONAL_RIVETLDADD@ @CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@;;
#  --python-libs)
#    echo @PYLIBDIR@;;
#  --load-mpi4py)
#    echo @LOADMPIFORPY@;;
#  --incdir)
#    echo ${SHERPA_INCLUDE_PATH:-@INCLUDEDIR@};;
#  --libdir)
#    echo ${SHERPA_LIBRARY_PATH:-@LIBDIR@};;
#  --datadir)
#    echo ${SHERPA_SHARE_PATH:-@DATADIR@};;
#  --ldflags)
#    echo -L${SHERPA_LIBRARY_PATH:-@LIBDIR@} -lSherpaMain 
#@CONDITIONAL_HEPMC2LIBS@ @CONDITIONAL_HEPMC3LIBS@ @CONDITIONAL_RIVETLDADD@ @CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@;;
#  --cppflags)
#    echo -I${SHERPA_INCLUDE_PATH:-@INCLUDEDIR@} 
#@CONDITIONAL_HEPMC2INCS@ @CONDITIONAL_HEPMC3INCS@ @CONDITIONAL_RIVETCPPFLAGS@ @CONDITIONAL_ROOTFLAGS@ @CONDITIONAL_LHAPDFINCS@;;
#  --cxx)

file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS_PACKAGES.H
     CONTENT 
"/* ATOOLS/Org/CXXFLAGS_PACKAGES.H.  Generated from CXXFLAGS_PACKAGES.H.in by configure.  */\n\
#ifndef CXXFLAGS_PACKAGES_H\n\
#define CXXFLAGS_PACKAGES_H\n\
\n\
/* Using HEPMC3 */\n\
#define USING__HEPMC3 \"1\"\n\
\n\
/* HepMC3 with ROOT support */\n\
/* #undef USING__HEPMC3__ROOT */\n\
\n\
/* Using HEPMC2 */\n\
#define USING__HEPMC2 \"1\"\n\
\n\
/* HepMC::IO_GenEvent available */\n\
#define USING__HEPMC2__IOGENEVENT \"1\"\n\
\n\
/* HepMC::Units available */\n\
#define USING__HEPMC2__UNITS \"1\"\n\
\n\
/* HepMCDefs.h available */\n\
#define USING__HEPMC2__DEFS \"1\"\n\
\n\
/* using gzip */\n\
#define USING__GZIP \"1\"\n\
\n\
/* using sqlite */\n\
/* #undef USING__SQLITE */\n\
\n\
/* hztool found */\n\
#define USING__HZTOOL \"1\"\n\
\n\
/* cernlib found */\n\
#define USING__CERNLIB \"1\"\n\
\n\
/* using LHAPDF */\n\
#define LHAPDF_PATH \"/usr\"\n\
\n\
/* using ROOT */\n\
#define USING__ROOT \"1\"\n\
\n\
/* using Rivet2 */\n\
/* #undef USING__RIVET2 */\n\
\n\
/* using Rivet3 */\n\
#define USING__RIVET3 \"1\"\n\
\n\
/* using BLACKHAT */\n\
#define USING__BLACKHAT \"1\"\n\
\n\
/* BLACKHAT directory */\n\
#define BLACKHAT_PATH \"/usr\"\n\
\n\
/* OpenLoops installation prefix */\n\
#define OPENLOOPS_PREFIX \"${OPENLOOPS_PREFIX}\"\n\
\n\
/* GoSam installation prefix */\n\
/* #undef GOSAM_PREFIX */\n\
\n\
/* MadLoop installation prefix */\n\
/* #undef MADLOOP_PREFIX */\n\
\n\
/* Recola installation prefix */\n\
#define RECOLA_PREFIX \"/usr\"\n\
\n\
/* using MCFM */\n\
/* #undef USING__MCFM */\n\
\n\
/* MCFM directory */\n\
/* #undef MCFM_PATH */\n\
\n\
/* using Delphes */\n\
/* #undef USING__DELPHES */\n\
\n\
/* using EWSud */\n\
/* #undef USING__EWSud */\n\
\n\
#endif"
     @ONLY)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/LundTools/pythia-6.4.18.f.in ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/LundTools/pythia-6.4.18.f @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/makelibs.in ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/makelibs @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/make2scons.in ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/make2scons @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/Sherpa-config.in ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/Run/Sherpa-config @ONLY)



file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS.H
     CONTENT 
"/* ATOOLS/Org/CXXFLAGS.H.  Generated from CXXFLAGS.H.in by configure.  */\n\
/* Architecture identified as Darwin MacOS */\n\
/* #undef ARCH_DARWIN */\n\
\n\
/* Architecture identified as Linux */\n\
#define ARCH_LINUX \"1\"\n\
\n\
/* Architecture identified as Unix */\n\
/* #undef ARCH_UNIX */\n\
\n\
/* Define to dummy `main' function (if any) required to link to the Fortran\n\
   libraries. */\n\
/* #undef FC_DUMMY_MAIN */\n\
\n\
/* Define if F77 and FC dummy `main' functions are identical. */\n\
/* #undef FC_DUMMY_MAIN_EQ_F77 */\n\
\n\
/* MPI support */\n\
/* #undef USING__MPI */\n\
\n\
/* Thread support */\n\
/* #undef USING__Threading */\n\
\n\
/* ld path name */\n\
#define LD_PATH_NAME \"${LD_LIBRARY_PATH}\"\n\
\n\
/* library suffix */\n\
#define LIB_SUFFIX \"${CMAKE_SHARED_LIBRARY_SUFFIX}\"\n\
\n\
/* binreloc activation */\n\
#define ENABLE_BINRELOC \"1\"\n\
\n\
/* Sherpa installation prefix */\n\
#define SHERPA_PREFIX \"${CMAKE_INSTALL_DIR}\"\n\
\n\
/* Sherpa include directory */\n\
#define SHERPA_INCLUDE_PATH \"${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC\"\n\
\n\
/* Sherpa library directory */\n\
#define SHERPA_LIBRARY_PATH \"${CMAKE_INSTALL_LIBDIR}/SHERPA-MC/\"\n\
\n\
/* Sherpa data directory */\n\
#define SHERPA_SHARE_PATH \"${CMAKE_INSTALL_DATADIR}/SHERPA-MC/\"\n\
\n\
/* Sherpa subversion */\n\
#define SHERPA_SUBVERSION \"0.0alpha1\"\n\
\n\
/* Sherpa version */\n\
#define SHERPA_VERSION \"3\"\n\
\n\
/* Sherpa version name */\n\
#define SHERPA_NAME \"Dhaulagiri\"\n\
\n\
/* using Pythia */\n\
#define USING__PYTHIA \"1\"\n\
\n\
/* HEPEVT common block size */\n\
#define HEPEVT_CB_SIZE ${HEPEVT_CB_SIZE}\n\
\n\
/* Using colour */\n\
#define USING__COLOUR \"1\"\n\
\n\
#define PYTHON_LIBS \"${SHERPA-MC_Python_SITEARCH}\"\n"
     @ONLY)
install(FILES 
      ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS.H 
      ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS_PACKAGES.H DESTINATION 
      ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC/ATOOLS/Org/  COMPONENT devel )
