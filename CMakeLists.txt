cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)
###########################Package@#####################################
project(SHERPA-MC LANGUAGES C  CXX Fortran )
#The complicated way to set the version is used to assure the simple 
#replacement of version variable with sed will not break the CMakeLists.txt
if (NOT __CL_VERSION__)
  set(EXTERNAL_VERSION_STRING "2.2.13")
else()
  set(EXTERNAL_VERSION_STRING "__CL_VERSION__.1")
endif()
string(REPLACE "."  ";" EXTERNAL_VERSION_LIST "${EXTERNAL_VERSION_STRING}")
list(GET EXTERNAL_VERSION_LIST 0 SHERPA_VERSION_MAJOR)
list(GET EXTERNAL_VERSION_LIST 1 SHERPA_VERSION_MINOR)
list(GET EXTERNAL_VERSION_LIST 2 SHERPA_VERSION_PATCH)
set(PROJECT_VERSION ${SHERPA_VERSION_MAJOR}.${SHERPA_VERSION_MINOR}.${SHERPA_VERSION_PATCH})
SET(SHERPA_VERSION ${SHERPA_VERSION_MAJOR}.${SHERPA_VERSION_MINOR})
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
option(SHERPA_ENABLE_PYTHON                "Enables building of tests." ON)
option(SHERPA_ENABLE_TEST                "Enables building of tests." ON)
option(SHERPA_ENABLE_MEMORYTEST          "Adds memory check tests." OFF)
option(SHERPA_POSITION_INDEPENDENT_CODE  "Position independent code for static libraries." OFF)
option(SHERPA_BUILD_SHARED               "Build shared libraries." ON)
option(SHERPA_BUILD_STATIC               "Build static libraries." OFF)
option(PYTHIA_SUPPORT               "Pythia support" ON)
option(GZIP_SUPPORT                "GZIP support" ON)
###########################Options handeling#####################################
if(SHERPA_ENABLE_TEST)
  ENABLE_TESTING()
  if (SHERPA_ENABLE_MEMORYTEST)
    find_program(VG NAMES valgrind REQUIRED)
    set(VGOPTIONS "--tool=memcheck --leak-check=full")
  endif()
endif()
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


find_package(HepMC3)
message(STATUS "SHERPA-MC: HEPMC3_VERSION=${HEPMC3_VERSION} HEPMC3_LIBRARIES=${HEPMC3_LIBRARIES}  HEPMC3_INCLUDE_DIR=${HEPMC3_INCLUDE_DIR}") 
find_package(HepMC2)
message(STATUS "SHERPA-MC: HEPMC2_VERSION=${HEPMC2_VERSION} HEPMC2_LIBRARIES=${HEPMC2_LIBRARIES}  HEPMC2_INCLUDE_DIR=${HEPMC2_INCLUDE_DIR}") 
find_package(LHAPDF REQUIRED)
message(STATUS "SHERPA-MC: LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS}")





add_subdirectory(AHADIC++)
add_subdirectory(AMEGIC++)
add_subdirectory(AMISIC++)
add_subdirectory(ATOOLS)

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
add_subdirectory(PDF/LHAPDF)  
add_subdirectory(PDF/Main )
add_subdirectory(PDF/NNPDF)  
add_subdirectory(PDF/SAL ) 
add_subdirectory(PDF/SASG)
add_subdirectory(RECONNECTIONS)
add_subdirectory(REMNANTS)
add_subdirectory(SHRiMPS)


add_subdirectory(SHERPA)


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
              EXCLUDE  PATTERN "*\.c" EXCLUDE )

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
                            # PATTERN "*Cluster_Algorithm.H" EXCLUDE 
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
#                         PATTERN Sherpa.tree.yaml EXCLUDE 
                         PATTERN Examples/Makefile EXCLUDE 
                         PATTERN Makefile.am EXCLUDE 
                         PATTERN Makefile.in EXCLUDE 
                         PATTERN Makefile.in EXCLUDE 
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
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/makelibs   ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/MadLoop/sconsloops DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PHASIC++/Main/plot_stats.sh  ${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/make2scons   
 ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/Sherpa-config
 DESTINATION  ${CMAKE_INSTALL_BINDIR} )


set(Sherpa_SRCS_C Sherpa.C Filter.C)
list(TRANSFORM Sherpa_SRCS_C PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Main/")
add_executable(Sherpa  ${Sherpa_SRCS_C} ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA//Run/Main.C)
target_link_libraries(Sherpa PRIVATE ATOOLS METOOLS ${LHAPDF_LIBRARIES} LHAPDF ${HEPMC3_LIBRARIES} ${HEPMC2_LIBRARIES} SHERPA PHASIC++ REMNANTS HADRONS++ PHOTONS++
AMISIC++ PHASIC++ BEAM SHRiMPS Main EXTRA_XS RECONNECTIONS  zip MODEL COMIX
)
target_include_directories(Sherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
install(TARGETS Sherpa DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libs)
#target_include_directories(Sherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/packlib/kuip)
#set_target_properties(Sherpa PROPERTIES POSITION_INDEPENDENT_CODE ${CERNLIB_POSITION_INDEPENDENT_CODE} OUTPUT_NAME Sherpa${COMPSUFFIX})




if (SHERPA_ENABLE_PYTHON)
add_subdirectory(AddOns/Python)

endif()























include(CMakePackageConfigHelpers)

set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/SHERPA/cmake)

#configure_package_config_file(cmake/Templates/SHERPAConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake/SHERPAConfig.cmake
#                              INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}
#                              PATH_VARS CMAKE_INSTALL_INCLUDEDIR
#                              CMAKE_INSTALL_LIBDIR)
            
#write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/cmake/SHERPAConfig-version.cmake COMPATIBILITY SameMajorVersion VERSION ${PROJECT_VERSION})


#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/SHERPAConfig.cmake
#              ${CMAKE_CURRENT_BINARY_DIR}/cmake/SHERPAConfig-version.cmake
#              DESTINATION ${CONFIG_INSTALL_DIR} COMPONENT devel)
              
