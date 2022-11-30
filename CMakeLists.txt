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

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE PATH "Default prefix" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
message(STATUS "SHERPA: CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
###########################Options#####################################
option(SHERPA_ENABLE_ANALYSIS                "Enables building of Analyses." OFF)
option(SHERPA_ENABLE_BINRELOC                "Enables binreloc." OFF)
option(SHERPA_ENABLE_BLACKHAT                "Enables building with BlackHat." OFF)
option(SHERPA_ENABLE_DIHIGGS                 "Enables building of DIHIGGS." OFF)
option(SHERPA_ENABLE_EWSUD                   "Enables building of EWSUD." OFF)
option(SHERPA_ENABLE_EXAMPLES                "Enables installation of EXAMPLES." ON)
option(SHERPA_ENABLE_GOSAM                   "Enables building with GoSam." OFF)
option(SHERPA_ENABLE_GZIP                    "Enables building with GZIP." OFF)
option(SHERPA_ENABLE_HEPMC2                  "Enables building with HepMC2." OFF)
option(SHERPA_ENABLE_HEPMC3                  "Enables building of HepMC3." OFF)
option(SHERPA_ENABLE_HEPMC3_ROOT             "Enables building of HepMC3+ROOT Interface ." OFF)
option(SHERPA_ENABLE_HZTOOL                  "Enables building with HZTool." OFF)
option(SHERPA_ENABLE_INTERNAL_PDFS           "Enables internal PDFs." OFF)
option(SHERPA_ENABLE_INSTALL_LIBZIP          "Enables internal LibZip." OFF)
option(SHERPA_ENABLE_INSTALL_LHAPDF          "Enables internal LHAPDF." OFF)
option(SHERPA_ENABLE_LHOLE                   "Enables building of LH OLE." OFF)
option(SHERPA_ENABLE_MADLOOP                 "Enables building with MadLoop." OFF)
option(SHERPA_ENABLE_MANUAL                  "Enables building of manul." OFF)
option(SHERPA_ENABLE_MCFM                    "Enables building with MCFM." OFF)
option(SHERPA_ENABLE_MPI                     "Enables building with MPI." OFF)
option(SHERPA_ENABLE_OPENLOOPS               "Enables building with OpenLoops." OFF)
option(SHERPA_ENABLE_PYTHIA6                 "Enables building with Pythia6 (embedded)." OFF)
option(SHERPA_ENABLE_PYTHIA8                 "Enables building with Pythia8." OFF)
option(SHERPA_ENABLE_PYTHON                  "Enables building of Python interface." OFF)
option(SHERPA_ENABLE_RECOLA                  "Enables building with Recola." OFF)
option(SHERPA_ENABLE_RIVET                   "Enables building with Rivet." OFF)
option(SHERPA_ENABLE_ROOT                    "Enables building with ROOT." OFF)
option(SHERPA_ENABLE_THREADING               "Enables building with threading." OFF)
option(SHERPA_ENABLE_TESTING                 "Enables testing." OFF)
option(SHERPA_ENABLE_UFO                     "Enables installation of UFO interfaces." OFF)
option(SHERPA_ENABLE_VERSIONING              "Enables versioned installation." OFF)
###########################Options handeling############################
set( allenables ANALYSIS BINRELOC BLACKHAT DIHIGGS EWSUD EXAMPLES GOSAM GZIP HEPMC2 HEPMC3 HEPMC3_ROOT HZTOOL
INTERNAL_PDFS INSTALL_LIBZIP INSTALL_LHAPDF 
 LHOLE MADLOOP MANUAL MCFM MPI OPENLOOPS PYTHIA6 PYTHIA8 PYTHON RECOLA RIVET ROOT THREADING TESTING UFO VERSIONING)
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
message(STATUS "SHERPA: CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}")
if (SHERPA_ENABLE_VERSIONING)
  set (VERSIONING "-${EXTERNAL_VERSION_STRING}")
else()
  set (VERSIONING "")
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
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING}")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING}")
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
#-flat_namespace 
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}  -single_module -undefined dynamic_lookup")
  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}  -single_module -undefined dynamic_lookup")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -single_module -undefined dynamic_lookup")
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

find_program(GIT NAMES git)
if (GIT)
  execute_process(COMMAND ${GIT} rev-parse --abbrev-ref HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITURL 
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GIT} rev-parse  HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITREV 
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "SHERPA: GIT=${GIT} GITREV=${GITREV} GITURL=${GITURL}")
else()
  message(STATUS "SHERPA: GIT IS NOT AVAILABLE!")
  set(GITURL "unknownurl")
  set(GITREV "unknownrevision")
endif()


macro(sherpa_create_git_info name)
  if (GIT)
    execute_process(COMMAND ${GIT} status -s --untracked-files=no .
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE GITCLEAN 
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
   set(GITCLEAN "")
  endif()
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
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/${f}"  temp)
    string(APPEND content "${temp}")
  endforeach()
  string(REPLACE "${PROJECT_SOURCE_DIR}/" "" GITTAG "${CMAKE_CURRENT_SOURCE_DIR}")
  string(MD5 MDF "${content}")
  set(newgitinfo "#include \"ATOOLS/Org/Git_Info.H\"\nstatic ATOOLS::Git_Info initializer\n(\"${GITTAG}\",\"${GITURL}\",\"${GITREV}${GITREVSUFFIX}\",\"${MDF}\");\n") 
  set(oldgitinfo "")
  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C oldgitinfo)
  endif()
  if ( NOT "${newgitinfo}" STREQUAL "${oldgitinfo}") 
    file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C "${newgitinfo}") 
  endif()
endmacro(sherpa_create_git_info )

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

macro(sherpa_handle_shared_library  mylib mylibname)
  if(SHERPA_ENABLE_MPI)
     target_link_libraries(${mylib} PRIVATE MPI::MPI_CXX)
     if (SHERPA_ENABLE_FORTRAN)
       target_link_libraries(${mylib} PRIVATE MPI::MPI_Fortran)
     endif()
  endif()
  target_include_directories(${mylib}  PRIVATE ${FREETYPE_INCLUDE_DIRS})
  install(TARGETS ${mylib} DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING} COMPONENT libs) 
  set_target_properties(${mylib}        PROPERTIES POSITION_INDEPENDENT_CODE ON  SOVERSION 0.0.0)
  set_target_properties(${mylib}        PROPERTIES DEFINE_SYMBOL "")
  set_target_properties(${mylib}        PROPERTIES OUTPUT_NAME ${mylibname})
  target_compile_features(${mylib} PRIVATE cxx_std_11)
endmacro(sherpa_handle_shared_library  mylib mylibname)

macro(sherpa_handle_executable  mylib )
  if(SHERPA_ENABLE_MPI)
     target_link_libraries(${mylib} PRIVATE MPI::MPI_CXX)
     if (SHERPA_ENABLE_FORTRAN)
       target_link_libraries(${mylib} PRIVATE MPI::MPI_Fortran)
     endif()
  endif()
  set_target_properties(${mylib} PROPERTIES  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>
                                             LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>
                                             RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR}/$<0:>
                                               )
   if(NOT (CMAKE_SYSTEM_NAME MATCHES Darwin))                                         
     target_link_options(${mylib} PRIVATE "LINKER:-no-as-needed")                                          
   endif()
   target_compile_features(${mylib} PRIVATE cxx_std_11)
endmacro(sherpa_handle_executable  mylib )


if (NOT SHERPA_ENABLE_INSTALL_LIBZIP)
  find_package(LibZip)
  if (NOT LibZip_FOUND)
    message(FATAL_ERROR "SHERPA: LibZip not found! Terminating!\n To use installed LibZip add to the cmake configuration -DLibZip_ROOT_DIR=/installation/path/to/libzip or use -DSHERPA_ENABLE_INSTALL_LIBZIP=ON to download and compile LibZip")
  endif()
else()
  message(STATUS "SHERPA: LibZip will be downloaded and compiled")
  include(ExternalProject)
  include(FetchContent)
  ExternalProject_Add(downloadedlibzip
    GIT_REPOSITORY https://github.com/nih-at/libzip
    GIT_TAG v1.9.2
    STEP_TARGETS build
    GIT_SHALLOW 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    PREFIX ${PROJECT_BINARY_DIR}/EXTERNALSRC
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/EXTERNALINSTALL/libzip -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DBUILD_DOC=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TOOLS=OFF -DBUILD_REGRESS=OFF
    BUILD_BYPRODUCTS ${PROJECT_BINARY_DIR}/EXTERNALINSTALL/libzip/lib/libzip${CMAKE_SHARED_LIBRARY_SUFFIX}
  )
  unset(LibZip_ROOT)
  set(LibZip_DIR ${PROJECT_BINARY_DIR}/EXTERNALINSTALL/libzip)
  set(LibZip_VERSION 1.9.2)
  set(LibZip_LIBRARIES ${LibZip_DIR}/lib/libzip${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(LibZip_INCLUDE_DIRS ${LibZip_DIR}/include)
  file(MAKE_DIRECTORY ${LibZip_DIR}/include)
  add_library(LibZip::LibZip UNKNOWN IMPORTED)
  set_target_properties(LibZip::LibZip PROPERTIES
        IMPORTED_LOCATION "${LibZip_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${LibZip_INCLUDE_DIRS}"
    )
  if (CMAKE_SYSTEM_NAME MATCHES Darwin)
    install(FILES ${LibZip_DIR}/lib/libzip.5.5${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
    install(FILES ${LibZip_DIR}/lib/libzip.5${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
  else()
    install(FILES ${LibZip_LIBRARIES}.5.5 DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
    install(FILES ${LibZip_LIBRARIES}.5 DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)  
  endif()
  install(FILES ${LibZip_LIBRARIES} DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
  install(FILES ${LibZip_INCLUDE_DIRS}/zip.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
  add_dependencies(LibZip::LibZip PRIVATE downloadedlibzip)
endif()
message(STATUS "SHERPA: LibZip_VERSION=${LibZip_VERSION} LibZip_LIBRARIES=${LibZip_LIBRARIES} LibZip_INCLUDE_DIR=${LibZip_INCLUDE_DIR}")

if (NOT SHERPA_ENABLE_INSTALL_LHAPDF)
  find_package(LHAPDF)
  if (NOT LHAPDF_FOUND)
    message(FATAL_ERROR "SHERPA: LHAPDF not found! Terminating!\n To use installed LHAPDF add to the cmake configuration -DLHAPDF_ROOT_DIR=/installation/path/to/lhapdf or use -DSHERPA_ENABLE_INSTALL_LHAPDF=ON to download and compile LHAPDF")
  endif()
else()
  message(STATUS "SHERPA: LHAPDF will be downloaded and compiled")
  include(ExternalProject)
  include(FetchContent)
  ExternalProject_Add(downloadedlhapdf
    URL https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.3.tar.gz
    UPDATE_COMMAND ""
    #LHAPDF is  bugged
    #One should not have doxygen in the path!
    #PATCH_COMMAND  "sed -i 's/ac_cv_path_DOXYGEN/DOXYGEN/g' configure"
    BUILD_IN_SOURCE     1
    PREFIX ${PROJECT_BINARY_DIR}/EXTERNALSRC
    CONFIGURE_COMMAND ./configure --disable-python --without-doxygen --prefix=${PROJECT_BINARY_DIR}/EXTERNALINSTALL/lhapdf
    BUILD_COMMAND make -j 8
    INSTALL_COMMAND make install
    BUILD_BYPRODUCTS ${PROJECT_BINARY_DIR}/EXTERNALINSTALL/lhapdf/lib/libLHAPDF${CMAKE_SHARED_LIBRARY_SUFFIX}
  )
  unset(LHAPDF_ROOT)
  set(LHAPDF_DIR ${PROJECT_BINARY_DIR}/EXTERNALINSTALL/lhapdf)
  set(LHAPDF_VERSION 6.5.3)
  set(LHAPDF_LIBRARIES ${LHAPDF_DIR}/lib/libLHAPDF${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(LHAPDF_INCLUDE_DIRS ${LHAPDF_DIR}/include)
  file(MAKE_DIRECTORY ${LHAPDF_DIR}/include)
  set(LHAPDF_DIR ${PROJECT_BINARY_DIR}/EXTERNALINSTALL/lhapdf)
  add_library(LHAPDF::LHAPDF UNKNOWN IMPORTED)
  set_target_properties(LHAPDF::LHAPDF PROPERTIES
        IMPORTED_LOCATION "${LHAPDF_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${LHAPDF_INCLUDE_DIRS}"
    )
  install(FILES ${LHAPDF_LIBRARIES} DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
  install(DIRECTORY ${LHAPDF_INCLUDE_DIRS}/LHAPDF DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
  install(DIRECTORY ${LHAPDF_DIR}/share/LHAPDF  DESTINATION ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT libs)
  add_dependencies(LHAPDF::LHAPDF PRIVATE downloadedlhapdf)
  
endif()
message(STATUS "SHERPA: LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS}")

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
add_subdirectory(PDF/LHAPDF)
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
  if (SHERPA_ENABLE_HEPMC2)
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

set( IDIRS  AHADIC++ AMEGIC++ AMISIC++ ATOOLS BEAM COMIX CSSHOWER++ DIM DIRE EXTAMP EXTRA_XS HADRONS++
METOOLS MCATNLO MODEL PHASIC++ PHOTONS++ PDF RECONNECTIONS REMNANTS SHRiMPS SHERPA)

if (SHERPA_ENABLE_MANUAL)
  add_subdirectory(Manual)
endif()

if (SHERPA_ENABLE_PYTHON)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/Python 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/AddOns/ COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE 
                             )
endif()

if (SHERPA_ENABLE_ANALYSIS)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/Analysis 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/AddOns/ COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" 
                             PATTERN "*/Scripts" EXCLUDE 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE 
                             )
endif()

if (SHERPA_ENABLE_EWSUD)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/EWSud
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/AddOns/ COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" 
                             EXCLUDE  PATTERN "*\.c" EXCLUDE 
                             EXCLUDE  PATTERN "*/Tests" EXCLUDE 
                             )
endif()
if (SHERPA_ENABLE_PYTHON OR SHERPA_ENABLE_DIHIGGS OR SHERPA_ENABLE_UFO)
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
#This directory is unconditional install, but requires Python for extra files to be installed
add_subdirectory(MODEL/UFO)




foreach ( IDIR ${IDIRS} )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${IDIR} 
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}  COMPONENT devel  
              FILES_MATCHING PATTERN "*.h" PATTERN "*.H" PATTERN "*.icc"  
                             PATTERN "*\.c" EXCLUDE 
                             PATTERN "*/fjcore.hh" 
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
                             PATTERN "*makefile*" EXCLUDE  
                             )    
endforeach (  IDIR ${IDIRS}  )
#Workaround for the Mac OS
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ATOOLS/Org/Node.C  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/ATOOLS/Org/)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ATOOLS/Org/Getter_Function.C  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/ATOOLS/Org/)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/LundTools/Dummy.C      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/LundTools/)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/LundTools/Lund_Interface.C      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/LundTools/)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ATOOLS/Math/Cluster_Algorithm.C      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC${VERSIONING}/ATOOLS/Math/)
if (SHERPA_ENABLE_EXAMPLES)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Examples DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} 
                         COMPONENT doc 
                         PATTERN test.py.in EXCLUDE 
                         PATTERN CI EXCLUDE 
                         PATTERN Examples/Makefile EXCLUDE 
                         PATTERN Makefile.am EXCLUDE 
                         PATTERN Makefile.in EXCLUDE 
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
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CT14/CT14Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc PATTERN "*ct14ll*" EXCLUDE  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/GRS/GRSGrid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/SAL/SALGrid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )
endif()
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/HEPEVT DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Math/Sobol/21201${SOBOL_EXT} DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING}/Sobol COMPONENT doc  )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/COMIX/Main/init_nlo.sh DESTINATION ${CMAKE_INSTALL_BINDIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/get_bibtex DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/plot_graphs DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2 DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )
add_custom_target( unZip ALL)
#TAR IS INTERNAL CMAKE COMMAND!
add_custom_command(TARGET unZip PRE_BUILD
       COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/NNPDF31_nnlo_as_0118_mc
       COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2
       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2
       COMMENT "UnpackingNNPDF31_nnlo_as_0118_mc.tar.bz2"
       VERBATIM)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/NNPDF31_nnlo_as_0118_mc DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc )

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/HADRONS++/Decaydata.zip DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc  )
install(FILES    ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/MadLoop/sconsloops DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} COMPONENT doc PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PHASIC++/Main/plot_stats.sh  
                ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/Run/Sherpa${VERSIONING}-config
                DESTINATION  ${CMAKE_INSTALL_BINDIR}  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

set(Sherpa_SRCS_C Sherpa.C Filter.C)
list(TRANSFORM Sherpa_SRCS_C PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Main/")
add_executable(Sherpa  ${Sherpa_SRCS_C} ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA//Run/Main.C)
target_link_libraries(Sherpa PRIVATE LibZip::LibZip)
sherpa_handle_executable(Sherpa)
set_target_properties(Sherpa        PROPERTIES OUTPUT_NAME Sherpa${VERSIONING})
target_link_libraries(Sherpa PRIVATE    
SherpaMain SherpaInitialization SherpaSingleEvents SherpaSoftPhysics SherpaPerturbativePhysics LundTools SherpaTools 
ShrimpsMain ShrimpsEvents ShrimpsLadders ShrimpsBeamRemnants ShrimpsXsecs ShrimpsEikonals ShrimpsTools 
HadronsMain HadronsMEs HadronsCurrents HadronsPSs 
Reconnections 
AhadicMain AhadicTools AhadicFormation AhadicDecays 
Amisic AmisicTools AmisicPerturbative 
PhotonsMain PhotonsTools PhotonsPhaseSpace PhotonsMEs 
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

if (SHERPA_ENABLE_RIVET)
  find_package(Rivet REQUIRED)
  message(STATUS "SHERPA: RIVET_VERSION=${RIVET_VERSION} RIVET_FOUND=${RIVET_FOUND}  RIVET_LIBRARIES=${RIVET_LIBRARIES} \
  RIVET_INCLUDE_DIRS=${RIVET_INCLUDE_DIRS} RIVET_DATA_PATH=${RIVET_DATA_PATH} RIVET_ANALYSIS_PATH=${RIVET_ANALYSIS_PATH} \
  Rivet_HEPMC2_FOUND=${Rivet_HEPMC2_FOUND}  Rivet_HEPMC3_FOUND=${Rivet_HEPMC3_FOUND} \
  RIVET_CONFIG_LIBS_STRING=${RIVET_CONFIG_LIBS_STRING}  RIVET_CONFIG_LIBS=${RIVET_CONFIG_LIBS}  RIVET_CONFIG_LIB_DIRS=${RIVET_CONFIG_LIB_DIRS} \
  RIVET_CONFIG_CPPFLAGS_STRING=${RIVET_CONFIG_CPPFLAGS_STRING}  RIVET_CONFIG_CPPFLAGS_DIRS=${RIVET_CONFIG_CPPFLAGS_DIRS} ")
  add_subdirectory(AddOns/Rivet)
endif()

if (SHERPA_ENABLE_PYTHIA8)
  find_package(Pythia8 8.220 REQUIRED)
  message(STATUS "SHERPA: PYTHIA8_VERSION=${PYTHIA8_VERSION} PYTHIA8_LIBRARIES=${PYTHIA8_LIBRARIES} PYTHIA8_INCLUDE_DIRS=${PYTHIA8_INCLUDE_DIRS}")
  add_subdirectory(AddOns/Pythia)
endif()

if (SHERPA_ENABLE_ROOT)
  #This is a hack to prevent looking for the nlohmannjson
  set(ROOT_builtin_nlohmannjson_FOUND OK)
  find_package(ROOT REQUIRED) 
  set(CMAKE_CXX_FLAGS_BEFORE_ROOT "${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS_BEFORE_ROOT "${CMAKE_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS_BEFORE_ROOT "${CMAKE_Fortran_FLAGS}")
  include(${ROOT_USE_FILE})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BEFORE_ROOT}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_BEFORE_ROOT}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_BEFORE_ROOT}")
  find_program(ROOT_CONFIG_EXE NAMES root-config PATHS ${ROOT_BINDIR}  NO_DEFAULT_PATH )
  
  set(ROOT_CONFIG_CPPFLAGS_STRING)
  if (ROOT_CONFIG_EXE)
    execute_process(COMMAND ${ROOT_CONFIG_EXE} --cflags
                  OUTPUT_VARIABLE ROOT_CONFIG_CPPFLAGS_STRING
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  message(STATUS "SHERPA: ROOT_VERSION=${ROOT_VERSION} ROOT_CONFIG_EXE=${ROOT_CONFIG_EXE} ROOT_LIBRARIES=${ROOT_LIBRARIES} ROOT_INCLUDE_DIRS=${ROOT_INCLUDE_DIRS} ROOT_CXX_FLAGS=${ROOT_CXX_FLAGS}")
  add_subdirectory(AddOns/Root)
endif()

if (SHERPA_ENABLE_BLACKHAT)
  find_package(BlackHat REQUIRED)
  message(STATUS "SHERPA: BLACKHAT_VERSION=${BLACKHAT_VERSION} BLACKHAT_PATH=${BLACKHAT_PATH}  BLACKHAT_LIBRARIES=${BLACKHAT_LIBRARIES} BLACKHAT_INCLUDE_DIRS=${BLACKHAT_INCLUDE_DIRS}")
  add_subdirectory(AddOns/BlackHat)
endif()

if (SHERPA_ENABLE_MADLOOP)
  if (MADLOOP_DIR)
    find_path(MADLOOP_PREFIX  mg5_aMC PATHS ${MADLOOP_DIR} PATH_SUFFIXES bin  REQUIRED NO_DEFAULT_PATH)
  else()
    find_program( MG5AMC NAMES  mg5_aMC REQUIRED)
    if (MG5AMC)
      get_filename_component (MADLOOP_PREFIX_BIN ${MG5AMC} DIRECTORY)
      get_filename_component (MADLOOP_PREFIX ${MADLOOP_PREFIX_BIN} DIRECTORY)
    endif()
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
  message(STATUS "SHERPA: RECOLA_VERSION=${RECOLA_VERSION} RECOLA_PREFIX=${RECOLA_PREFIX} RECOLA_INCLUDE_DIR=${RECOLA_INCLUDE_DIR}")
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
get_filename_component (CXX "${CMAKE_CXX_COMPILER}" NAME)
#Here we just use the typical flags used for the compilation
get_target_property(CXXSTANDARD Sherpa CXX_STANDARD)
if (NOT CXXSTANDARD)
  set(CXXSTANDARD "")
else()
  set(CXXSTANDARD "-std=c++${CXXSTANDARD}")
endif()
set(CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${CXXSTANDARD}")
message(STATUS "Using for the script CXX=${CXX} CXXFLAGS=${CXXFLAGS}")
set (CXXCPP "${CXX} -E")
set (FCFLAGS ${CMAKE_Fortran_FLAGS})
set (AM_CXXFLAGS ${CMAKE_CXX_FLAGS})
get_filename_component (CC "${CMAKE_C_COMPILER}" NAME)
set (CFLAGS ${CMAKE_C_FLAGS})
if (SHERPA_ENABLE_FORTRAN)
  get_filename_component (FC "${CMAKE_Fortran_COMPILER}" NAME)
endif()
set (INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC)
set (LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC)
set (BINDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
set (DATADIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/SHERPA-MC)
set (prefix ${CMAKE_INSTALL_PREFIX})
set (PACKAGE_VERSION ${SHERPA_VERSION})
set (host "${CMAKE_HOST_SYSTEM_PROCESSOR}-pc-somehost")
if (NOT HEPEVT_CB_SIZE)
  set(HEPEVT_CB_SIZE  4000)
endif()
set( LD_LIBRARY_PATH "LD_LIBRARY_PATH")
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
  set( LD_LIBRARY_PATH "DYLD_LIBRARY_PATH")
endif()
set (PYLIBDIR ${SHERPA_Python_SITEARCH})

if (SHERPA_ENABLE_ROOT)
  set(CONDITIONAL_ROOTFLAGS "${ROOT_CONFIG_CPPFLAGS_STRING}")
endif()
if (SHERPA_ENABLE_HEPMC2)
  set(CONDITIONAL_HEPMC2LIBS "-L${HEPMC2_LIBRARY_DIR} -lHepMC")
  set(CONDITIONAL_HEPMC2INCS "-I${HEPMC2_INCLUDE_DIRS}")
else()
  set(CONDITIONAL_HEPMC2LIBS )
  set(CONDITIONAL_HEPMC2INCS )
endif()
if (SHERPA_ENABLE_HEPMC3)
  set(CONDITIONAL_HEPMC3INCS "-I${HEPMC3_INCLUDE_DIR}")
  get_filename_component (HEPMC3_LIBRARY_DIR ${HEPMC3_LIB} DIRECTORY)
  if (SHERPA_ENABLE_HEPMC3_ROOT)
    get_filename_component (HEPMC3_ROOTIO_LIBRARY_DIR ${HEPMC3_ROOTIO_LIB} DIRECTORY)
    set(CONDITIONAL_HEPMC3LIBS "-L${HEPMC3_LIBRARY_DIR} -lHepMC3 -L${HEPMC3_ROOTIO_LIBRARY_DIR} -lHepMC3rootIO")
  else()
    set(CONDITIONAL_HEPMC3LIBS "-L${HEPMC3_LIBRARY_DIR} -lHepMC3")
  endif()
else()
  set(CONDITIONAL_HEPMC3INCS )
  set(CONDITIONAL_HEPMC3LIBS )
endif()

set(CONDITIONAL_LHAPDFLIBS  "-L${LHAPDF_LIBRARY_DIR} -lLHAPDF")
set(CONDITIONAL_LHAPDFINCS  "${LHAPDF_CONFIG_CPPFLAGS_STRING}")

if (SHERPA_ENABLE_RIVET)
  set(CONDITIONAL_RIVETLDADD " ${RIVET_CONFIG_LIBS_STRING}")
  set(CONDITIONAL_RIVETCPPFLAGS "${RIVET_CONFIG_CPPFLAGS_STRING}")
else()
  set(CONDITIONAL_RIVETLDADD )
  set(CONDITIONAL_RIVETCPPFLAGS )
endif()
########################################################################
if (SHERPA_ENABLE_RIVET AND RIVET_VERSION VERSION_LESS 3.0.0)
  set(RIVET2 1)
endif()
if (SHERPA_ENABLE_RIVET AND RIVET_VERSION VERSION_GREATER 2.99.99)
  set(RIVET3 1)
endif()
set(CXXFLAGS_PACKAGES_H_CONTENT
"/* ATOOLS/Org/CXXFLAGS_PACKAGES.H.  Generated from CXXFLAGS_PACKAGES.H.in by configure.  */\n\
#ifndef CXXFLAGS_PACKAGES_H\n\
#define CXXFLAGS_PACKAGES_H\n\n\
/* Using HEPMC3 */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_HEPMC3}>,#define USING__HEPMC3 \"1\",/* #undef USING__HEPMC3 */>\n\n\
/* HepMC3 with ROOT support */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_HEPMC3_ROOT}>,#define USING__HEPMC3__ROOT \"1\",/* #undef USING__HEPMC3__ROOT */>\n\n\
/* Using HEPMC2 */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_HEPMC2}>,#define USING__HEPMC2 \"1\",/* #undef USING__HEPMC2 */>\n\n\
/* HepMC::IO_GenEvent available */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_HEPMC2}>,#define USING__HEPMC2__IOGENEVENT \"1\",/* #undef USING__HEPMC2__IOGENEVENT */>\n\n\
/* HepMC::Units available */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_HEPMC2}>,#define USING__HEPMC2__UNITS \"1\",/* #undef USING__HEPMC2__UNITS */>\n\n\
/* HepMCDefs.h available */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_HEPMC2}>,#define USING__HEPMC2__DEFS \"1\",/* #undef USING__HEPMC2__DEFS */>\n\n\
/* using gzip */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_GZIP}>,#define USING__GZIP \"1\",/* #undef USING__GZIP */>\n\n\
/* hztool found */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_HZTOOL}>,#define USING__HZTOOL \"1\",/* #undef USING__HZTOOL */>\n\n\
/* cernlib found */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_CERNLIB}>,#define USING__CERNLIB \"1\",/* #undef USING__CERNLIB */>\n\n\
/* using LHAPDF */\n\
#define LHAPDF_PATH \"${LHAPDF_PATH}\"\n\n\
/* using ROOT */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_ROOT}>,#define USING__ROOT \"1\",/* #undef USING__ROOT */>\n\n\
/* using Rivet2 */\n\
$<IF:$<BOOL:${RIVET2}>,#define USING__RIVET2 \"1\",/* #undef USING__RIVET2 */>\n\n\
/* using Rivet3 */\n\
$<IF:$<BOOL:${RIVET3}>,#define USING__RIVET3 \"1\",/* #undef USING__RIVET3 */>\n\n\
/* using BLACKHAT */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_BLACKHAT}>,#define USING__BLACKHAT \"1\",/* #undef USING__BLACKHAT */>\n\n\
/* BLACKHAT directory */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_BLACKHAT}>,#define BLACKHAT_PATH \"${BLACKHAT_PATH}\",/* #undef BLACKHAT_PATH */>\n\n\
/* OpenLoops installation prefix */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_OPENLOOPS}>,#define OPENLOOPS_PREFIX \"${OPENLOOPS_PREFIX}\",/* #undef OPENLOOPS_PREFIX */>\n\n\
/* GoSam installation prefix */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_GOSAM}>,#define GOSAM_PREFIX \"${GOSAM_PREFIX}\",/* #undef GOSAM_PREFIX */>\n\n\
/* MadLoop installation prefix */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_MADLOOP}>,#define MADLOOP_PREFIX \"${MADLOOP_PREFIX}\",/* #undef MADLOOP_PREFIX */>\n\n\
/* Recola installation prefix */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_RECOLA}>,#define RECOLA_PREFIX \"${RECOL_PREFIX}\",/* #undef RECOLA_PREFIX */>\n\n\
/* using MCFM */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_MCFM}>,#define USING__MCFM \"1\",/* #undef USING__MCFM */>\n\n\
/* MCFM directory */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_MCFM}>,#define MCFM_PATH \"${MCFM_PATH}\",/* #undef MCFM_PATH */>\n\n\
/* using EWSud */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_EWSUD}>,#define USING__EWSud \"1\",/* #undef USING__EWSud */>\n\n\
#endif")
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS_PACKAGES.H CONTENT "${CXXFLAGS_PACKAGES_H_CONTENT}")
########################################################################
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/LundTools/pythia-6.4.18.f.in ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/LundTools/pythia-6.4.18.f @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/makelibs.in ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/makelibs @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/makelibs    DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc PERMISSIONS  OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/make2scons.in ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/make2scons @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/Sherpa-config.in ${CMAKE_CURRENT_BINARY_DIR}/SHERPA/Run/Sherpa-config @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AMEGIC++/Main/make2scons   DESTINATION  ${CMAKE_INSTALL_BINDIR} PERMISSIONS  OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

if (SHERPA_ENABLE_UFO)
  set (PYTHON "${Python_EXECUTABLE}")
  set (PYLIBDIR "${SHERPA_Python_SITEARCH}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/MODEL/UFO/sconstruct_template.in ${CMAKE_CURRENT_BINARY_DIR}/MODEL/UFO/sconstruct_template @ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/MODEL/UFO/Sherpa-generate-model.in ${CMAKE_CURRENT_BINARY_DIR}/MODEL/UFO/Sherpa-generate-model @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MODEL/UFO/Sherpa-generate-model DESTINATION  ${CMAKE_INSTALL_BINDIR} PERMISSIONS  OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MODEL/UFO/sconstruct_template DESTINATION ${SHERPA_Python_SITEARCH}/ufo_interface/  COMPONENT python PERMISSIONS  OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
endif()
########################################################################
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(DARWIN_SYS 1)
endif()
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
  set(LINUX_SYS 1)
endif()
if (UNIX AND ( NOT  CMAKE_SYSTEM_NAME MATCHES "Linux")  AND (NOT (CMAKE_SYSTEM_NAME MATCHES Darwin)))
  set(UNIX_SYS 1)
endif()
set(CXXFLAGS_H_CONTENT "/* ATOOLS/Org/CXXFLAGS.H.  Generated from CXXFLAGS.H.in by configure.  */\n\
/* Architecture identified as Darwin MacOS */\n\
$<IF:$<BOOL:${DARWIN_SYS}>,#define ARCH_DARWIN \"1\",/* #undef ARCH_DARWIN */>\n\n\
/* Architecture identified as Linux */\n\
$<IF:$<BOOL:${LINUX_SYS}>,#define ARCH_LINUX 1,/* #undef ARCH_LINUX */>\n\n\
/* Architecture identified as Unix */\n\
$<IF:$<BOOL:${UNIX_SYS}>,#define ARCH_UNIX 1,/* #undef ARCH_UNIX */>\n\n\
/* Define to dummy `main' function (if any) required to link to the Fortran\n   libraries. */\n\
/* #undef FC_DUMMY_MAIN */\n\n\
/* Define if F77 and FC dummy `main' functions are identical. */\n\
/* #undef FC_DUMMY_MAIN_EQ_F77 */\n\n\
/* MPI support */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_MPI}>,#define USING__MPI \"1\",/* #undef USING__MPI */>\n\n\
/* Thread support */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_THREADING}>,#define USING__Threading \"1\",/* #undef USING__Threading */>\n\n\
/* ld path name */\n\
#define LD_PATH_NAME \"${LD_LIBRARY_PATH}\"\n\n\
/* library suffix */\n\
#define LIB_SUFFIX \"${CMAKE_SHARED_LIBRARY_SUFFIX}\"\n\n\
/* binreloc activation */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_BINRELOC}>,#define ENABLE_BINRELOC \"1\",/* #undef ENABLE_BINRELOC */>\n\n\
/* Sherpa installation prefix */\n\
#define SHERPA_PREFIX \"${CMAKE_INSTALL_PREFIX}\"\n\n\
/* Sherpa include directory */\n\
#define SHERPA_INCLUDE_PATH \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC\"\n\n\
/* Sherpa library directory */\n\
#define SHERPA_LIBRARY_PATH \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/SHERPA-MC/\"\n\n\
/* Sherpa data directory */\n\
#define SHERPA_SHARE_PATH \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/SHERPA-MC/\"\n\n\
/* Sherpa subversion */\n\
#define SHERPA_SUBVERSION \"${SHERPA_VERSION_MINOR}.${SHERPA_VERSION_PATCH}\"\n\n\
/* Sherpa version */\n\
#define SHERPA_VERSION \"${SHERPA_VERSION_MAJOR}\"\n\n\
/* Sherpa version name */\n\
#define SHERPA_NAME \"Dhaulagiri\"\n\n\
/* using Pythia */\n\
$<IF:$<BOOL:${SHERPA_ENABLE_PYTHIA6}>,#define USING__PYTHIA \"1\",/* #undef USING__PYTHIA */>\n\n\
/* HEPEVT common block size */\n\
#define HEPEVT_CB_SIZE ${HEPEVT_CB_SIZE} \n\n\
/* Using colour */\n\
#define USING__COLOUR \"1\"\n\n\
#define PYTHON_LIBS \"${SHERPA_Python_SITEARCH}\"\n")
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ATOOLS/Org/CXXFLAGS.H CONTENT "${CXXFLAGS_H_CONTENT}")
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
