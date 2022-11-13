cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)
###########################Package@#####################################
project(SHERPA-MC LANGUAGES C  CXX)
#The complicated way to set the version is used to assure the simple 
#replacement of version variable with sed will not break the CMakeLists.txt
set(EXTERNAL_VERSION_STRING "3.0.0alpha1")
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
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH} )

if (APPLE)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
###########################Options#####################################
option(SHERPA_ENABLE_ANALYSIS                "Enables building of Analyses." ON)
option(SHERPA_ENABLE_BINRELOC                "Enables binreloc." ON)
option(SHERPA_ENABLE_BLACKHAT                "Enables building with BlackHat." ON)
option(SHERPA_ENABLE_DIHIGGS                 "Enables building of DIHIGGS." ON)
option(SHERPA_ENABLE_EWSUD                   "Enables building of EWSUD." ON)
option(SHERPA_ENABLE_EXAMPLES                "Enables installation of EXAMPLES." ON)
option(SHERPA_ENABLE_GOSAM                   "Enables building with GoSam." ON)
option(SHERPA_ENABLE_GZIP                    "Enables building with GZIP." ON)
option(SHERPA_ENABLE_HEPMC2                  "Enables building with HepMC2." ON)
option(SHERPA_ENABLE_HEPMC3                  "Enables building of HepMC3." ON)
option(SHERPA_ENABLE_HEPMC3_ROOT             "Enables building of HepMC3+ROOT Interface ." ON)
option(SHERPA_ENABLE_HZTOOL                  "Enables building with HZTool." ON)
option(SHERPA_ENABLE_INTERNAL_PDFS           "Enables internal PDFs." ON)
option(SHERPA_ENABLE_LHAPDF                  "Enables building With LHAPDF." ON)
option(SHERPA_ENABLE_LHOLE                   "Enables building of LH OLE." ON)
option(SHERPA_ENABLE_MADLOOP                 "Enables building with MadLoop." ON)
option(SHERPA_ENABLE_MANUAL                  "Enables building of manul." ON)
option(SHERPA_ENABLE_MCFM                    "Enables building with MCFM." ON)
option(SHERPA_ENABLE_MPI                     "Enables building with MPI." OFF)
option(SHERPA_ENABLE_OPENLOOPS               "Enables building with OpenLoops." ON)
option(SHERPA_ENABLE_PYTHIA6                 "Enables building with Pythia6 (embedded)." ON)
option(SHERPA_ENABLE_PYTHIA8                 "Enables building with Pythia8." ON)
option(SHERPA_ENABLE_PYTHON                  "Enables building of Python interface." ON)
option(SHERPA_ENABLE_RECOLA                  "Enables building with Recola." ON)
option(SHERPA_ENABLE_RIVET                   "Enables building with Rivet." ON)
option(SHERPA_ENABLE_ROOT                    "Enables building with ROOT." ON)
option(SHERPA_ENABLE_THREADING               "Enables building SHERPA_ENABLE_MANUALwith threading." ON)#BUGGED?
option(SHERPA_ENABLE_TESTING                 "Enables testing." ON)
###########################Options handeling############################
set( allenables ANALYSIS BINRELOC BLACKHAT DIHIGGS EWSUD GOSAM GZIP HEPMC2 HEPMC3 HEPMC3_ROOT HZTOOL 
INTERNAL_PDFS
LHAPDF LHOLE MADLOOP MANUAL MCFM MPI OPENLOOPS PYTHIA6 PYTHIA8 PYTHON RECOLA RIVET ROOT THREADING TESTING)
foreach (en IN LISTS allenables)
     message(STATUS "SHERPA: SHERPA_ENABLE_${en}=${SHERPA_ENABLE_${en}}")
endforeach()
if (SHERPA_ENABLE_PYTHIA6 OR SHERPA_ENABLE_INTERNAL_PDFS)
 set(SHERPA_ENABLE_FORTRAN ON)
 enable_language(Fortran)
 message(STATUS "SHERPA: SHERPA_ENABLE_FORTRAN=${SHERPA_ENABLE_FORTRAN}")
else()
 message(STATUS "SHERPA: SHERPA_ENABLE_FORTRAN=${SHERPA_ENABLE_FORTRAN}")
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
if (SHERPA_ENABLE_FORTRAN)
  get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  message(STATUS "SHERPA: Fortran_COMPILER_NAME=${Fortran_COMPILER_NAME}  CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}  CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
endif()
get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
message(STATUS "SHERPA: C_COMPILER_NAME=${C_COMPILER_NAME}  CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}   CMAKE_C_COMPILER_VERSION=${CMAKE_C_COMPILER_VERSION}")
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
message(STATUS "SHERPA: CXX_COMPILER_NAME=${CXX_COMPILER_NAME}  CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}   CMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}")

if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
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
message(STATUS "SHERPA: CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
set(CXX_CXX_FLAGS_TO_CHECK -Wno-tautological-constant-compare)
set ( NEW_CMAKE_CXX_FLAGS " ")
string(REPLACE " " ";" LIST_CMAKE_CXX_FLAGS "${CXX_CXX_FLAGS_TO_CHECK} ${CMAKE_CXX_FLAGS}")
foreach(fx IN LISTS LIST_CMAKE_CXX_FLAGS)
  if (fx)
    check_cxx_compiler_flag(${fx} CXX_COMPILER_SUPPORTS_${fx})
    if(CXX_COMPILER_SUPPORTS_${fx})
      list( APPEND NEW_CMAKE_CXX_FLAGS ${fx})
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${fx}" ) #neede for dependent flags
    endif()
  endif()
endforeach()
list(JOIN NEW_CMAKE_CXX_FLAGS " " CMAKE_CXX_FLAGS_X )
if (NOT CMAKE_FORCE_FLAGS)
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_X})
endif()
message(STATUS "SHERPA: CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
########################################################################

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
if (SHERPA_ENABLE_FORTRAN)
  include(CheckFortranCompilerFlag)
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
    set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_X})
  endif()
  message(STATUS "SHERPA: CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
endif()
########################################################################
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
add_definitions(-DHAVE_CONFIG_H)
add_definitions(-DPIC)

find_program(GIT NAMES git  REQUIRED)
execute_process(COMMAND ${GIT} rev-parse --abbrev-ref HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITURL 
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${GIT} rev-parse  HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITREV 
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "SHERPA: GIT=${GIT} GITREV=${GITREV} GITURL=${GITURL}")
macro(creategitinfo name b)
  execute_process(COMMAND ${GIT} status -s --untracked-files=no .
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITCLEAN 
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if ("${GITCLEAN}" STREQUAL "")
    set(GITREVSUFFIX "")
  else()
    set(GITREVSUFFIX "-dirty")
#Message only changes
    message(STATUS "SHERPA ${CMAKE_CURRENT_SOURCE_DIR}: GITCLEAN=${GITCLEAN} GITREVSUFFIX=${GITREVSUFFIX}")
  endif()
  set(allsources )
  foreach (f IN LISTS ${name}_SOURCES  ${name}_HEADERS )
    list(APPEND allsources ${f})
  endforeach()
  list (FILTER allsources EXCLUDE REGEX ".*CXXFLAGS.*")
  set(content "")
  foreach (f IN LISTS allsources )
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/${b}/${f}"  temp)
    string(APPEND content "${temp}")
  endforeach()
  string(MD5 MDF "${content}")
  set(newgitinfo "#include \"ATOOLS/Org/Git_Info.H\"\nstatic ATOOLS::Git_Info initializer\n(\"${GITTAG}\",\"${GITURL}\",\"${GITREV}${GITREVSUFFIX}\",\"${MDF}\");\n") 
  set(oldgitinfo "")
  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C oldgitinfo)
  endif()
  if ( NOT "${newgitinfo}" STREQUAL "${oldgitinfo}") 
    file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C "${newgitinfo}") 
  endif()
endmacro(creategitinfo )

if(SHERPA_ENABLE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "SHERPA: MPI_HOME=${MPI_HOME}")
  if (SHERPA_ENABLE_FORTRAN)
    message(STATUS "SHERPA: MPI_CXX_COMPILE_OPTIONS=${MPI_CXX_COMPILE_OPTIONS} MPI_C_COMPILE_OPTIONS=${MPI_C_COMPILE_OPTIONS} MPI_Fortran_COMPILE_OPTIONS=${MPI_Fortran_COMPILE_OPTIONS}")
    message(STATUS "SHERPA: MPI_CXX_COMPILE_DEFINITIONS=${MPI_CXX_COMPILE_DEFINITIONS} MPI_C_COMPILE_DEFINITIONS=${MPI_C_COMPILE_DEFINITIONS} MPI_Fortran_COMPILE_DEFINITIONS=${MPI_Fortran_COMPILE_DEFINITIONS}")
    message(STATUS "SHERPA: MPI_CXX_INCLUDE_DIRS=${MPI_CXX_INCLUDE_DIRS} MPI_C_INCLUDE_DIRS=${MPI_C_INCLUDE_DIRS} MPI_Fortran_INCLUDE_DIRS=${MPI_Fortran_INCLUDE_DIRS}")
    message(STATUS "SHERPA: MPI_CXX_LINK_FLAGS=${MPI_CXX_LINK_FLAGS} MPI_C_LINK_FLAGS=${MPI_C_LINK_FLAGS} MPI_Fortran_LINK_FLAGS=${MPI_Fortran_LINK_FLAGS}")
    include_directories(${MPI_Fortran_INCLUDE_DIRS})
  else()
    message(STATUS "SHERPA: MPI_CXX_COMPILE_OPTIONS=${MPI_CXX_COMPILE_OPTIONS} MPI_C_COMPILE_OPTIONS=${MPI_C_COMPILE_OPTIONS}")
    message(STATUS "SHERPA: MPI_CXX_COMPILE_DEFINITIONS=${MPI_CXX_COMPILE_DEFINITIONS} MPI_C_COMPILE_DEFINITIONS=${MPI_C_COMPILE_DEFINITIONS}")
    message(STATUS "SHERPA: MPI_CXX_INCLUDE_DIRS=${MPI_CXX_INCLUDE_DIRS} MPI_C_INCLUDE_DIRS=${MPI_C_INCLUDE_DIRS}")
    message(STATUS "SHERPA: MPI_CXX_LINK_FLAGS=${MPI_CXX_LINK_FLAGS} MPI_C_LINK_FLAGS=${MPI_C_LINK_FLAGS}")
  endif()
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  set (LOADMPIFORPY "from mpi4py import MPI")
else()
  set (LOADMPIFORPY "")
endif()

macro(sherpa_mpi_link_libraries  mylib)
  if(SHERPA_ENABLE_MPI)
     target_link_libraries(${mylib} PRIVATE MPI::MPI_CXX)
     if (SHERPA_ENABLE_FORTRAN)
       target_link_libraries(${mylib} PRIVATE MPI::MPI_Fortran)
     endif()
  endif()
endmacro(sherpa_mpi_link_libraries  mylib)

#One should have also zlib?
find_package(LibZip REQUIRED)
#FIND_LIBRARY(LIBZIP_LIBRARY NAMES zip PATHS "$ENV{LIB_DIR}/lib" "$ENV{LIB}" /usr/local/lib /usr/lib /usr/lib64 )
message(STATUS "SHERPA: LibZip_VERSION=${LibZip_VERSION} LibZip_LIBRARIES=${LibZip_LIBRARIES} LibZip_INCLUDE_DIR=${LibZip_INCLUDE_DIR}")

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
add_subdirectory(PDF/Electron) 
if (SHERPA_ENABLE_INTERNAL_PDFS) 
  add_subdirectory(PDF/GRS)
  add_subdirectory(PDF/GRV)
  add_subdirectory(PDF/CT14)
  add_subdirectory(PDF/SASG)
endif()
IF (SHERPA_ENABLE_LHAPDF)
  find_package(LHAPDF REQUIRED)
  message(STATUS "SHERPA: LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS}")
  add_subdirectory(PDF/LHAPDF)
endif()
add_subdirectory(ATOOLS)
add_subdirectory(PDF/Main)
add_subdirectory(PDF/NNPDF)  
add_subdirectory(PDF/SAL) 
add_subdirectory(RECONNECTIONS)
add_subdirectory(REMNANTS)
add_subdirectory(SHRiMPS)
IF (SHERPA_ENABLE_MCFM)
  find_package(MCFM REQUIRED)
  message(STATUS "SHERPA: MCFM_VERSION=${MCFM_VERSION} MCFM_LIBRARY=${MCFM_LIBRARY} MCFM_INCLUDE_DIR=${MCFM_INCLUDE_DIR}")
  add_subdirectory(AddOns/MCFM)
endif()

if (SHERPA_ENABLE_HEPMC3 OR SHERPA_ENABLE_HEPMC2)
  if (SHERPA_ENABLE_HEPMC3)
     find_package(HepMC2 REQUIRED)
     message(STATUS "SHERPA: HEPMC2_VERSION=${HEPMC2_VERSION} HEPMC2_LIBRARIES=${HEPMC2_LIBRARIES} HEPMC2_FIO_LIBRARIES=${HEPMC2_FIO_LIBRARIES} HEPMC2_INCLUDE_DIRS=${HEPMC2_INCLUDE_DIRS}")
  endif()
  if (SHERPA_ENABLE_HEPMC3)
    find_package(HepMC3 REQUIRED)
    message(STATUS "SHERPA: HEPMC3_VERSION=${HEPMC3_VERSION} HEPMC3_LIBRARIES=${HEPMC3_LIBRARIES}  HEPMC3_INCLUDE_DIR=${HEPMC3_INCLUDE_DIR}") 
    if (SHERPA_ENABLE_HEPMC3_ROOT)
      if (NOT HEPMC3_ROOTIO_LIB)
        message(FATAL_ERROR "SHERPA: HEPMC3_ROOTIO_LIB NOT FOUND")
      else()
        message(STATUS "SHERPA: HEPMC3_ROOTIO_LIB=${HEPMC3_ROOTIO_LIB}")
      endif()
    endif()
  endif()
  add_subdirectory(AddOns/HepMC)
endif()
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
SHERPA)

if (SHERPA_ENABLE_MANUAL)
  add_subdirectory(Manual)
endif()

if (SHERPA_ENABLE_PYTHON)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/Python 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC/AddOns/ COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" 
                             PATTERN "*.deps*" EXCLUDE   
                             PATTERN "*.lib*" EXCLUDE 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE 
                             )
endif()

if (SHERPA_ENABLE_ANALYSIS)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/Analysis 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC/AddOns/ COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" 
                             PATTERN "*.deps*" EXCLUDE   
                             PATTERN "*.lib*" EXCLUDE 
                             PATTERN "*/Scripts" EXCLUDE 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE 
                             )
endif()

if (SHERPA_ENABLE_EWSUD)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/EWSud
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC/AddOns/ COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" 
                             PATTERN "*.deps*" EXCLUDE   
                             PATTERN "*.lib*" EXCLUDE 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE 
                             EXCLUDE  PATTERN "*/Tests" EXCLUDE 
                             )
endif()

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
if (SHERPA_ENABLE_EXAMPLES)
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
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ATOOLS/Math/Sobol/21201 ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Math/Sobol/21201)

if (SHERPA_ENABLE_GZIP)
  find_program(GZIP NAMES gzip  REQUIRED)
  message(STATUS "SHERPA: GZIP=${GZIP}")
  set(SOBOL_EXT .gz)
  add_custom_target( ZipSobol ALL)
  add_custom_command(TARGET ZipSobol PRE_BUILD
       COMMAND ${GZIP} -f 21201
       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Math/Sobol
       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Math/Sobol/21201
       VERBATIM)
endif()
if (SHERPA_ENABLE_INTERNAL_PDFS)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CT14/CT14Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc PATTERN "*ct14ll*" EXCLUDE  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/GRS/GRSGrid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/SAL/SALGrid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
endif()
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/HEPEVT DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Math/Sobol/21201${SOBOL_EXT} DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC/Sobol COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/COMIX/Main/init_nlo.sh DESTINATION  ${CMAKE_INSTALL_BINDIR} )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/get_bibtex DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/plot_graphs DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2 DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
add_custom_target( unZip ALL)
#TAR IS INTERNAL CMAKE COMMAND!
add_custom_command(TARGET unZip PRE_BUILD
       COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/NNPDF31_nnlo_as_0118_mc
       COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2
       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2
       COMMENT "UnpackingNNPDF31_nnlo_as_0118_mc.tar.bz2"
       VERBATIM)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/NNPDF31_nnlo_as_0118_mc DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc )

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/HADRONS++/Decaydata.zip DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES    ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/MadLoop/sconsloops DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PHASIC++/Main/plot_stats.sh  
                ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/Run/Sherpa-config
                DESTINATION  ${CMAKE_INSTALL_BINDIR} )

set(Sherpa_SRCS_C Sherpa.C Filter.C)
list(TRANSFORM Sherpa_SRCS_C PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Main/")
add_executable(Sherpa  ${Sherpa_SRCS_C} ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA//Run/Main.C)
target_link_libraries(Sherpa PRIVATE LibZip::LibZip)

target_link_libraries(Sherpa PRIVATE    
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
sherpa_mpi_link_libraries(Sherpa)

if (SHERPA_ENABLE_HEPMC3)
  target_link_libraries(Sherpa PRIVATE SherpaHepMC3Output  ${HEPMC3_LIB} ${HEPMC3_SEARCH_LIB} )
  if (SHERPA_ENABLE_HEPMC3_ROOT)
    target_link_libraries(Sherpa PRIVATE SherpaHepMC3Output  ${HEPMC3_ROOTIO_LIB} )
  endif()
endif()
if (SHERPA_ENABLE_HEPMC2)
  target_link_libraries(Sherpa PRIVATE SherpaHepMCOutput  ${HEPMC2_LIBRARIES} )
endif()

target_include_directories(Sherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
install(TARGETS Sherpa DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libs)

add_subdirectory(AddOns/Weights)
add_subdirectory(AddOns/NNLO)
add_subdirectory(AddOns/Higgs)
add_subdirectory(AddOns/EWSud)
if (SHERPA_ENABLE_PYTHON OR SHERPA_ENABLE_DIHIGGS)
  if (NOT  SHERPA-MC_PYTHON_VERSION)
    find_package (Python 3 COMPONENTS  Development Interpreter)
  else()
    find_package (Python ${SHERPA-MC_PYTHON_VERSION} EXACT COMPONENTS  Development Interpreter)
  endif()
  if (NOT SHERPA_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
    set (SHERPA_Python_SITEARCH ${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/)
    message(STATUS "SHERPA: Python_VERSION_MAJOR.Python_VERSION_MINOR=${Python_VERSION_MAJOR}${Python_VERSION_MINOR} SHERPA_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR} defaults to ${SHERPA_Python_SITEARCH}")
  else()
    message(STATUS "SHERPA: Python_VERSION_MAJOR.Python_VERSION_MINOR=${Python_VERSION_MAJOR}${Python_VERSION_MINOR}Python_VERSION_MAJOR.Python_VERSION_MINOR=${Python_VERSION_MAJOR}${Python_VERSION_MINOR} SHERPA_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR} is set to to ${SHERPA_Python_SITEARCH}")
  endif()
  if (SHERPA_ENABLE_PYTHON)
    add_subdirectory(AddOns/Python)
  endif()
  if (SHERPA_ENABLE_DIHIGGS)
    add_subdirectory(AddOns/DiHiggsNLO)
  endif()
endif()

if (SHERPA_ENABLE_RIVET)
  find_package(Rivet REQUIRED)
  message(STATUS "SHERPA: RIVET_VERSION=${RIVET_VERSION} RIVET_FOUND=${RIVET_FOUND}  RIVET_LIBRARIES=${RIVET_LIBRARIES} \
  RIVET_INCLUDE_DIRS=${RIVET_INCLUDE_DIRS} RIVET_DATA_PATH=${RIVET_DATA_PATH} RIVET_ANALYSIS_PATH=${RIVET_ANALYSIS_PATH} \
  Rivet_HEPMC2_FOUND=${Rivet_HEPMC2_FOUND}  Rivet_HEPMC3_FOUND=${Rivet_HEPMC3_FOUND} \
  RIVET_CONFIG_LIBS_STRING=${RIVET_CONFIG_LIBS_STRING}  RIVET_CONFIG_LIBS=${RIVET_CONFIG_LIBS}  RIVET_CONFIG_LIB_DIRS=${RIVET_CONFIG_LIB_DIRS} \
  RIVET_CONFIG_CPPFLAGS_STRING=${RIVET_CONFIG_CPPFLAGS_STRING}  RIVET_CONFIG_CPPFLAGS_DIRS=${RIVET_CONFIG_CPPFLAGS_DIRS} \
  "
  )
  add_subdirectory(AddOns/Rivet)
endif()

if (SHERPA_ENABLE_PYTHIA8)
  find_package(Pythia8 8.220 REQUIRED)
  message(STATUS "SHERPA: PYTHIA8_VERSION=${PYTHIA8_VERSION} PYTHIA8_LIBRARIES=${PYTHIA8_LIBRARIES} PYTHIA8_INCLUDE_DIRS=${PYTHIA8_INCLUDE_DIRS}")
  add_subdirectory(AddOns/Pythia)
endif()

if (SHERPA_ENABLE_ROOT)
  add_subdirectory(AddOns/Root)
endif()

if (SHERPA_ENABLE_BLACKHAT)
  find_package(BlackHat REQUIRED)
  message(STATUS "SHERPA: BLACKHAT_VERSION=${BLACKHAT_VERSION} BLACKHAT_PATH=${BLACKHAT_PATH}  BLACKHAT_LIBRARIES=${BLACKHAT_LIBRARIES} BLACKHAT_INCLUDE_DIRS=${BLACKHAT_INCLUDE_DIRS}")
  add_subdirectory(AddOns/BlackHat)
endif()

if (SHERPA_ENABLE_MADLOOP)
  if (NOT MADLOOP_PREFIX)
    find_path(MADLOOP_PREFIX  mg5_aMC PATHS PATH_SUFFIXES bin  REQUIRED)
  else()
    find_program( MG5AMC NAMES  mg5_aMC PATHS ${MADLOOP_PREFIX} PATH_SUFFIXES bin  NO_DEFAULT_PATH REQUIRED)
  endif()
  message(STATUS "SHERPA: MADLOOP_PREFIX=${MADLOOP_PREFIX}")
  add_subdirectory(AddOns/MadLoop)
endif()


if (SHERPA_ENABLE_LHOLE)
  add_subdirectory(AddOns/LH_OLE)
endif()

if (SHERPA_ENABLE_OPENLOOPS)
  find_package(OpenLoops REQUIRED)
  message(STATUS "SHERPA: OPENLOOPS_PREFIX=${OPENLOOPS_PREFIX} OPENLOOPS_LIBRARY=${OPENLOOPS_LIBRARY}")
  add_subdirectory(AddOns/OpenLoops)
endif()

if (SHERPA_ENABLE_RECOLA)
  find_package(Recola REQUIRED)
  message(STATUS "SHERPA: RECOLA_VERSION=${RECOLA_VERSION} RECOLA_PREFIX=${RECOLA_PREFIX}")
  add_subdirectory(AddOns/Recola)
endif()

if (SHERPA_ENABLE_HZTOOL)
  find_package(HZTool REQUIRED)
  message(STATUS "SHERPA: HZTOOL_LIBRARY=${HZTOOL_LIBRARY} HZTOOL_INCLUDE_DIR=${HZTOOL_INCLUDE_DIR}")
  find_package(CERNLIB REQUIRED)
  message(STATUS "SHERPA: CERN_kernlib_DYNAMIC_LIBRARY=${CERN_kernlib_DYNAMIC_LIBRARY} CERN_mathlib_DYNAMIC_LIBRARY=${CERN_mathlib_DYNAMIC_LIBRARY} CERN_packlib_DYNAMIC_LIBRARY=${CERN_packlib_DYNAMIC_LIBRARY}")
  find_package(Motif REQUIRED)
  message(STATUS "SHERPA: MOTIF_INCLUDE_DIR=${MOTIF_INCLUDE_DIR} MOTIF_LIBRARIES=${MOTIF_LIBRARIES}")
  add_subdirectory(AddOns/HZTool)
endif()

if (SHERPA_ENABLE_GOSAM)
  find_package(GoSam REQUIRED)
  message(STATUS "SHERPA: GOSAM_VERSION=${GOSAM_VERSION} GOSAM_PATH=${GOSAM_PREFIX}")
  add_subdirectory(AddOns/GoSam)
endif()

IF (SHERPA_ENABLE_ANALYSIS)
  add_subdirectory(AddOns/Analysis)
endif()

########################################################################
# Config files
set (CXX ${CMAKE_CXX_COMPILER})
set (CXXFLAGS ${CMAKE_CXX_FLAGS})
set (CXXCPP "${CMAKE_CXX_COMPILER} -E")
set (FCFLAGS ${CMAKE_Fortran_FLAGS})
set (AM_CXXFLAGS ${CMAKE_CXX_FLAGS})
set (CC ${CMAKE_C_COMPILER})
set (CFLAGS ${CMAKE_C_FLAGS})
set (FC ${CMAKE_Fortran_COMPILER})
set (INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC)
set (LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC)
set (BINDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
set (DATADIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR})
set (prefix ${CMAKE_INSTALL_PREFIX})
set (PACKAGE_VERSION ${SHERPA_VERSION})
set (host "${CMAKE_HOST_SYSTEM_PROCESSOR}-pc-somehost")
if (NOT HEPEVT_CB_SIZE)
  set(HEPEVT_CB_SIZE  4000)
endif()
set( LD_LIBRARY_PATH "LD_LIBRARY_PATH")
if (Darwin)
  set( LD_LIBRARY_PATH "DYLD_LIBRARY_PATH")
endif()
set (PYLIBDIR ${SHERPA_Python_SITEARCH})
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
########################################################################
set(CXXFLAGS_PACKAGES_H_CONTENT
"/* ATOOLS/Org/CXXFLAGS_PACKAGES.H.  Generated from CXXFLAGS_PACKAGES.H.in by configure.  */\n\
#ifndef CXXFLAGS_PACKAGES_H\n\
#define CXXFLAGS_PACKAGES_H\n\
\n")
if (SHERPA_ENABLE_HEPMC3)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* Using HEPMC3 */\n#define USING__HEPMC3 \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* Using HEPMC3 */\n/* #undef USING__HEPMC3 */\n\n")
endif()
if (SHERPA_ENABLE_HEPMC3_ROOT)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMC3 with ROOT support */\n#define USING__HEPMC3__ROOT \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMC3 with ROOT support */\n/* #undef USING__HEPMC3__ROOT */\n\n")
endif()
if (SHERPA_ENABLE_HEPMC2)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* Using HEPMC2 */\n#define USING__HEPMC2 \"1\"\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMC::IO_GenEvent available */\n#define USING__HEPMC2__IOGENEVENT \"1\"\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMC::Units available */\n#define USING__HEPMC2__UNITS \"1\"\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMCDefs.h available */\n#define USING__HEPMC2__DEFS \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* Using HEPMC2 */\n/* #undef USING__HEPMC2 */\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMC::IO_GenEvent available */\n/* #undef USING__HEPMC2__IOGENEVENT */\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMC::Units available */\n/* #undef USING__HEPMC2__UNITS */\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* HepMCDefs.h available */\n/* #undef USING__HEPMC2__DEFS */\n\n")
endif()

if (SHERPA_ENABLE_GZIP)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using gzip */\n#define USING__GZIP \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using gzip */\n/* #undef USING__GZIP */\n\n")
endif()

if (SHERPA_ENABLE_SQLITE)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using sqlite */\n#define USING__SQLITE \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using sqlite */\n/* #undef USING__SQLITE */\n\n")
endif()

if (SHERPA_ENABLE_HZTOOL)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* hztool found */\n#define USING__HZTOOL \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* hztool found */\n/* #undef USING__HZTOOL */\n\n")
endif()

if (SHERPA_ENABLE_CERNLIB)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* cernlib found */\n#define USING__CERNLIB \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* cernlib found */\n/* #undef USING__CERNLIB */\n\n")
endif()

if (SHERPA_ENABLE_LHAPDF)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using LHAPDF */\n#define LHAPDF_PATH \"${LHAPDF_PATH}\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using LHAPDF */\n/* #undef LHAPDF_PATH */\n\n")
endif()

if (SHERPA_ENABLE_ROOT)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using ROOT */\n#define USING__ROOT \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using ROOT */\n/* #undef USING__ROOT */\n\n")
endif()

if (SHERPA_ENABLE_RIVET AND RIVET_VERSION VERSION_LESS 3.0.0)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using Rivet2 */\n#define USING__RIVET2 \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using Rivet2 */\n/* #undef USING__RIVET2 */\n\n")
endif()

if (SHERPA_ENABLE_RIVET AND RIVET_VERSION VERSION_GREATER 2.99.99)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using Rivet3 */\n#define USING__RIVET3 \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using Rivet3 */\n/* #undef USING__RIVET3 */\n\n")
endif()

if (SHERPA_ENABLE_BLACKHAT)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using BLACKHAT */\n#define USING__BLACKHAT \"1\"\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* BLACKHAT directory */\n#define BLACKHAT_PATH \"${BLACKHAT_PATH}\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using BLACKHAT */\n/* #undef USING__BLACKHAT */\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* BLACKHAT directory */\n/* #undef BLACKHAT_PATH */\n\n")
endif()

if (SHERPA_ENABLE_OPENLOOPS)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* OpenLoops installation prefix */\n#define OPENLOOPS_PREFIX \"${OPENLOOPS_PREFIX}\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* OpenLoops installation prefix */\n/* #undef OPENLOOPS_PREFIX */\n\n")
endif()

if (SHERPA_ENABLE_GOSAM)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* GoSam installation prefix */\n#define GOSAM_PREFIX \"${GOSAM_PREFIX}\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* GoSam installation prefix */\n/* #undef GOSAM_PREFIX */\n\n")
endif()

if (SHERPA_ENABLE_MADLOOP)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* MadLoop installation prefix */\n#define MADLOOP_PREFIX \"${MADLOOP_PREFIX}\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* MadLoop installation prefix */\n/* #undef MADLOOP_PREFIX */\n\n")
endif()

if (SHERPA_ENABLE_RECOLA)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* Recola installation prefix */\n#define RECOLA_PREFIX \"${RECOLA_PREFIX}\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* Recola installation prefix */\n/* #undef RECOLA_PREFIX */\n\n")
endif()

if (SHERPA_ENABLE_MCFM)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using MCFM */\n#define USING__MCFM \"1\"\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* MCFM directory */\n#define MCFM_PATH \"${MCFM_PATH}\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using MCFM */\n/* #undef USING__MCFM */\n\n")
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* MCFM directory */\n/* #undef MCFM_PATH */\n\n")
endif()

if (SHERPA_ENABLE_DELPHES)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using Delphes */\n#define USING__DELPHES \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using Delphes */\n/* #undef USING__DELPHES */\n\n")
endif()

if (SHERPA_ENABLE_EWSUD)
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using EWSud */\n#define USING__EWSud \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_PACKAGES_H_CONTENT "/* using EWSud */\n/* #undef USING__EWSud */\n\n")
endif()
string(APPEND CXXFLAGS_PACKAGES_H_CONTENT  "#endif")
file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS_PACKAGES.H CONTENT "${CXXFLAGS_PACKAGES_H_CONTENT}" @ONLY)
########################################################################

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/LundTools/pythia-6.4.18.f.in ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/LundTools/pythia-6.4.18.f @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/makelibs.in ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/makelibs @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/makelibs    DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/make2scons.in ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/make2scons @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/Sherpa-config.in ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/Run/Sherpa-config @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/make2scons   DESTINATION  ${CMAKE_INSTALL_BINDIR} )

########################################################################
set(CXXFLAGS_H_CONTENT "/* ATOOLS/Org/CXXFLAGS.H.  Generated from CXXFLAGS.H.in by configure.  */\n")

if (Darwin)
  string(APPEND CXXFLAGS_H_CONTENT  "/* Architecture identified as Darwin MacOS */\n#define ARCH_DARWIN \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_H_CONTENT  "/* Architecture identified as Darwin MacOS */\n/* #undef ARCH_DARWIN */\n\n")
endif()

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
  string(APPEND CXXFLAGS_H_CONTENT  "/* Architecture identified as Linux */\n#define ARCH_LINUX \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_H_CONTENT  "/* Architecture identified as Linux */\n/* #undef ARCH_LINUX */\n\n")
endif()

if (UNIX AND ( NOT  CMAKE_SYSTEM_NAME MATCHES "Linux")  AND (NOT Darwin) )
  string(APPEND CXXFLAGS_H_CONTENT  "/* Architecture identified as Unix */\n#define ARCH_UNIX 1\n\n")
else()
  string(APPEND CXXFLAGS_H_CONTENT  "/* Architecture identified as Unix */\n/* #undef ARCH_UNIX */\n\n")
endif()

#This is just obsolete
string(APPEND CXXFLAGS_H_CONTENT  "\
/* Define to dummy `main' function (if any) required to link to the Fortran\n\
   libraries. */\n\
/* #undef FC_DUMMY_MAIN */\n\
\n\
/* Define if F77 and FC dummy `main' functions are identical. */\n\
/* #undef FC_DUMMY_MAIN_EQ_F77 */\n\n")

if (SHERPA_ENABLE_MPI)
  string(APPEND CXXFLAGS_H_CONTENT  "/* MPI support */\n#define USING__MPI \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_H_CONTENT  "/* MPI support */\n/* #undef USING__MPI */\n\n")
endif()

if (SHERPA_ENABLE_THREADING)
  string(APPEND CXXFLAGS_H_CONTENT  "/* Thread support */\n#define USING__Threading \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_H_CONTENT  "/* Thread support */\n/* #undef USING__Threading */\n")
endif()

string(APPEND CXXFLAGS_H_CONTENT  "
/* ld path name */\n\
#define LD_PATH_NAME \"${LD_LIBRARY_PATH}\"\n\
\n\
/* library suffix */\n\
#define LIB_SUFFIX \"${CMAKE_SHARED_LIBRARY_SUFFIX}\"\n\
\n")

if (SHERPA_ENABLE_BINRELOC)
  string(APPEND CXXFLAGS_H_CONTENT "/* binreloc activation */\n#define ENABLE_BINRELOC \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_H_CONTENT "/* binreloc activation */\n/* #undef ENABLE_BINRELOC */\n\n")
endif()

string(APPEND CXXFLAGS_H_CONTENT  "\
/* Sherpa installation prefix */\n\
#define SHERPA_PREFIX \"${CMAKE_INSTALL_PREFIX}\"\n\
\n\
/* Sherpa include directory */\n\
#define SHERPA_INCLUDE_PATH \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC\"\n\
\n\
/* Sherpa library directory */\n\
#define SHERPA_LIBRARY_PATH \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC/\"\n\
\n\
/* Sherpa data directory */\n\
#define SHERPA_SHARE_PATH \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/SHERPA-MC/\"\n\
\n\
/* Sherpa subversion */\n\
#define SHERPA_SUBVERSION \"${SHERPA_VERSION_MINOR}.${SHERPA_VERSION_PATCH}\"\n\
\n\
/* Sherpa version */\n\
#define SHERPA_VERSION \"${SHERPA_VERSION_MAJOR}\"\n\
\n\
/* Sherpa version name */\n\
#define SHERPA_NAME \"Dhaulagiri\"\n\
\n")

if (SHERPA_ENABLE_PYTHIA6)
  string(APPEND CXXFLAGS_H_CONTENT  "/* using Pythia */\n#define USING__PYTHIA \"1\"\n\n")
else()
  string(APPEND CXXFLAGS_H_CONTENT  "/* using Pythia */\n/* #undef USING__PYTHIA */\n\n")
endif()

string(APPEND CXXFLAGS_H_CONTENT  "\
/* HEPEVT common block size */\n\
#define HEPEVT_CB_SIZE ${HEPEVT_CB_SIZE} \n\
\n\
/* Using colour */\n\
#define USING__COLOUR \"1\"\n\
\n\
#define PYTHON_LIBS \"${SHERPA_Python_SITEARCH}\"\n")

file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS.H CONTENT "${CXXFLAGS_H_CONTENT}" @ONLY)
########################################################################
install(FILES 
      ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS.H 
      ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS_PACKAGES.H DESTINATION 
      ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC/ATOOLS/Org/  COMPONENT devel )
########################################################################
if (SHERPA_ENABLE_TESTING)
  ENABLE_TESTING()
  if (SHERPA_ENABLE_EWSUD)
    add_subdirectory(AddOns/EWSud/Tests)
  endif()
endif()
