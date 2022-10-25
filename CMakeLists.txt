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
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
 set(LINUX 1)
endif()
if (CMAKE_SYSTEM_NAME MATCHES "SunOS")
 set(SunOS 1)
endif()
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
option(SHERPA_ENABLE_TEST                "Enables building of tests." ON)
option(SHERPA_ENABLE_MEMORYTEST          "Adds memory check tests." OFF)
option(SHERPA_POSITION_INDEPENDENT_CODE  "Position independent code for static libraries." OFF)
option(SHERPA_BUILD_SHARED               "Build shared libraries." ON)
option(SHERPA_BUILD_STATIC               "Build static libraries." OFF)
set(SHERPA_BUILD_STATIC ON)  #This is always ON, but it is convenient to have it like this from the point of view of logic.
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
########################################################################
if (APPLE)
  find_program ( SED NAMES gsed REQUIRED)
  find_program ( CPP NAMES clang REQUIRED)
endif()
if (LINUX)
  find_program ( SED NAMES sed REQUIRED)
  find_program ( CPP NAMES cpp REQUIRED)
endif()
if (SunOS)
  find_program ( SED NAMES gsed sed REQUIRED)
  find_program ( CPP NAMES cpp gcpp clang REQUIRED)
endif()
message(STATUS "SHERPA: SED=${SED}")
if (CMAKE_C_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(C_C_FLAGS_TO_CHECK       
                               "-fcommon"
                               "-fexceptions"
                               "-trigraphs"
                               "-mfpmath=387"    
                               "-Wno-cpp"  
                               "-Wno-conversion"  
                               "-Wno-uninitialized"  
                               "-Wno-implicit-int"  
                               "-Wno-int-to-pointer-cast"  
                               "-Wno-implicit-function-declaration" 
                               "-Wno-incompatible-pointer-types" 
                               "-Wno-return-type" 
                               "-Wno-trigraphs" 
                               "-Wno-unused-label" 
                               "-Wno-unused-function" 
                               "-Wno-unused-dummy-argument"
                               "-Wno-missing-braces"
                               "-Wno-switch"
                               "-Wno-format-overflow"
                               "-Wno-unused-result"
                               "-Wno-array-bounds"
                               "-Wno-restrict"
                               "-Wno-tautological-constant-compare"
                               "-Wno-constant-logical-operand"
                               "-O0"
                               "-g"
                               )
endif()
########################################################################
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(C_C_FLAGS_TO_CHECK       "-fcommon"
                               "--param=ssp-buffer-size=4"
                               "-fexceptions"
                               "-fasynchronous-unwind-tables"
                               "-fcf-protection"
                               "-trigraphs" 
                               "-mfpmath=387"                               
                               "-grecord-gcc-switches"
                               "-pipe"
                               "-Wall"
                               "-Wno-trigraphs" 
                               "-Wno-pointer-to-int-cast"
                               "-Wno-unused-but-set-variable" 
                               "-Wno-unused-function" 
                               "-Wno-missing-braces"
                               "-Wno-switch"
                               "-Wno-array-bounds"
                               "-Wno-conversion"  
                               "-Wno-uninitialized"  
                               "-Wno-int-to-pointer-cast"  
                               "-Wno-implicit-function-declaration" 
                               "-Wformat -Werror=format-security"
                               "-O0"
                               "-g"
                               )
if (CMAKE_C_COMPILER_VERSION VERSION_GREATER 4.0)
     list(APPEND C_C_FLAGS_TO_CHECK  
                                     "-fstack-protector-strong" 
                                     "-fstack-clash-protection"
                                     )
endif()
if (CMAKE_C_COMPILER_VERSION VERSION_GREATER 4.5)
    list(APPEND C_C_FLAGS_TO_CHECK 
                         "-Wno-restrict"
                         "-Wno-unused-result"
                         "-Wno-format-overflow"
                         "-Wno-unused-dummy-argument"
                         "-Wno-maybe-uninitialized"
                         "-Wno-stringop-truncation"
                         "-Wno-incompatible-pointer-types"
                         "-Wno-cpp"
                         "-fstack-protector"
                         )
endif()
endif()
if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(C_C_FLAGS_TO_CHECK "-Wno-cpp"
                         "-Wall"
                         "-fcommon"
                         "-diag-disable=117" # non-void function should return a value
                         "-diag-disable=140" # Too many arguments in function call
                         "-diag-disable=165" # Too few arguments in function call
                         "-diag-disable=167" # incompatible parameter
                         "-diag-disable=177"  # labels
                         "-diag-disable=180" # incompatible format
                         "-diag-disable=181"  # format
                         "-diag-disable=266" # implicit declaration
                         "-diag-disable=303" # explicit type is missing
                         "-diag-disable=556"  # value of type "long *" cannot be assigned to an entity of type "int *"
                         "-diag-disable=810" # conversion from char to int
                         "-diag-disable=1224" # _FORTIFY_SOURCE requires compiling with optimization (-O)
                         "-diag-disable=2193" # null argument 
                         "-O0"
                         "-mp1"
                         "-g"
                             )
endif()
if (CMAKE_C_COMPILER_ID STREQUAL "NVHPC")
  set(C_C_FLAGS_TO_CHECK       "-Wno-cpp"
                               "-fcommon"  
                               "-Wno-conversion"  
                               "-Wno-uninitialized"  
                               "-Wno-implicit-int"  
                               "-Wno-int-to-pointer-cast"  
                               "-Wno-implicit-function-declaration" 
                               "-Wno-incompatible-pointer-types" 
                               "-Wno-return-type" 
                               "-trigraphs" 
                               "-Wno-trigraphs" 
                               "-Wno-unused-label" 
                               "-Wno-unused-function" 
                               "-Wno-unused-dummy-argument"
                               "-Wno-missing-braces"
                               "-Wno-switch"
                               "-Wno-format-overflow"
                               "-Wno-unused-result"
                               "-Wno-array-bounds"
                               "-Wno-restrict"
                               "-flto=no"
                               "-O0" 
                               "-g"  
                               "-fstack-protector"
                               )
endif()
########################################################################
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(F_Fortran_FLAGS_TO_CHECK 
                               "--param=ssp-buffer-size=4"
                               "-fno-backslash" 
                               "-ffixed-line-length-72" 
                               "-std=legacy"
                               "-fcommon"
                               "-fno-automatic" 
                               "-fno-second-underscore" 
                               "-mfpmath=387" 
                               "-fno-align-commons" 
                               "-fno-aggressive-loop-optimizations" 
                               "-finit-local-zero"
                               "-fallow-argument-mismatch"
                               "-fno-range-check"
                               "-fno-bounds-check"
                               "-fbacktrace"
                               "-Wno-cpp" 
                               "-Wno-uninitialized" 
                               "-Wno-character-truncation" 
                               "-Wno-line-truncation" 
                               "-Wno-lto-type-mismatch"
                               "-Wno-integer-division"
                               "-Wno-intrinsic-shadow"
                               "-Wno-argument-mismatch"
                               "-Wno-aggressive-loop-optimizations"
                               "-Wno-unused-label" 
                               "-Wno-unused-function" 
                               "-Wno-unused-variable" 
                               "-Wno-maybe-uninitialized" 
                               "-Wno-unused-dummy-argument"
                               "-O0"
                               "-g"
                               "-fallow-invalid-boz" #needed for patchy4/p4lib/jartpx.F
                               )
#Not clear starting from which version we can rely on the stack protector
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 4.5)
    list(APPEND F_Fortran_FLAGS_TO_CHECK  "-fstack-protector" )
  endif()
  if (Fortran_COMPILER_NAME STREQUAL g77)
    SET(CMAKE_C_USE_RESPONSE_FILE_FOR_OBJECTS 0)
    SET(CMAKE_Fortran_USE_RESPONSE_FILE_FOR_OBJECTS 0)
  else()
    list(APPEND F_Fortran_FLAGS_TO_CHECK    
                               "-Wno-conversion" 
                               "-Wno-implicit-function-declaration"
                               "-Wno-return-type"
   )
  endif()
########################################################################
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(F_Fortran_FLAGS_TO_CHECK 
                              "-Wall"
                              "-fcommon"
                              "-g"
                              "-diag-disable6178" # return value is not defined
#                              "-diag-disable7615" 
                              "-diag-disable7713" # Statement function has not been used
                              "-diag-disable7784" # Symbol in BLOCK DATA program unit is not in a COMMON block.
                              "-diag-disable8291" # Recommended relationship between field width
                              "-diag-disable8577" # The scale factor (k) ..
                              "-O0"
                              "-traceback"
                              "-noalign" 
                              "-noauto" 
                              "-zero" 
                              "-nbs"
                              "-mp1"
   )
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(F_Fortran_FLAGS_TO_CHECK 
                               "-fcommon" 
                               "-traceback"
                               "-Mnobounds"
                               "-Mnostandard"
                               "-Mnobackslash" 
                               "-Mpreprocess" 
                               "-Mnosave" 
                               "-Mnofree"
                               "-Wno-cpp" 
                               "-O0" 
                               "-g"
   )
elseif()
  message(FATAL "SHERPA: Don't know how to set Fortran FLAGS for  ${Fortran_COMPILER_NAME} with ID ${Fortran_COMPILER_ID}")
endif ()
if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
endif()

#Linux has here -Wl,-z,relro -Wl,--as-needed  -Wl,-z,now -specs=/usr/lib/rpm/redhat/redhat-hardened-ld -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1
set(CMAKE_EXE_LINKER_FLAGS "")
set(CMAKE_MODULE_LINKER_FLAGS  "")
set(CMAKE_SHARED_LINKER_FLAGS  "")
set(CMAKE_STATIC_LINKER_FLAGS  "")

set(CMAKE_EXE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}  "")
set(CMAKE_MODULE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}  "")
set(CMAKE_SHARED_LINKER_FLAGS_${CMAKE_BUILD_TYPE}  "")
set(CMAKE_STATIC_LINKER_FLAGS_${CMAKE_BUILD_TYPE}  "")

set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}  "")
set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}  "")

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
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_X})
endif()
message(STATUS "SHERPA : CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
########################################################################
include_directories(${PROJECT_SOURCE_DIR})

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
#add_subdirectory(PDF)
add_subdirectory(RECONNECTIONS)
add_subdirectory(REMNANTS)
add_subdirectory(SHRiMPS)
add_subdirectory(SHERPA)

#find_package(ZIP)

  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Examples DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC 
  COMPONENT doc 
  FILES_MATCHING
                         PATTERN test.py EXCLUDE 
                         PATTERN CI EXCLUDE 
                         PATTERN Userhook EXCLUDE 
                         PATTERN HPC EXCLUDE 
                         PATTERN NLO_EW_* EXCLUDE 
                         PATTERN NLO_QCD_* EXCLUDE 
                         PATTERN RemnantTests EXCLUDE 
                         PATTERN VBF_Interference EXCLUDE 
                         PATTERN FixedOrder_NNLO EXCLUDE 
#                         PATTERN Sherpa.tree.yaml EXCLUDE 
                         PATTERN Makefile EXCLUDE 
                         PATTERN Makefile.am EXCLUDE 
                         PATTERN Makefile.in EXCLUDE 
                          )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/MRST/MRST01Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/MRST/MRST04Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/MRST/MRST99Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
 # install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118mc DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CT10/CT10Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CT12/CT12Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CT14/CT14Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CTEQ/CTEQ6Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/CTEQ/CTEQ66Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PDF/MSTW/MSTW08Grid DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AddOns/HEPEVT DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ATOOLS/Math/Sobol/21201.gz DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC/Sobol COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/COMIX/Main/init_nlo.sh DESTINATION  ${CMAKE_INSTALL_BINDIR} )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/get_bibtex DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/SHERPA/Run/plot_graphs DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PDF/NNPDF/NNPDF31_nnlo_as_0118_mc.tar.bz2 DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/AHADRONS++/Decaydata.zip DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/AMEGIC++/Main/makelibs DESTINATION  ${CMAKE_INSTALL_DATADIR}/SHERPA-MC COMPONENT doc  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PHASIC++/Main/plot_stats.sh
 DESTINATION  ${CMAKE_INSTALL_BINDIR} )
  install(FILES 
               ${CMAKE_CURRENT_SOURCE_DIR}/AUTHORS 
               ${CMAKE_CURRENT_SOURCE_DIR}/COPYING 
               ${CMAKE_CURRENT_SOURCE_DIR}/README 
               ${CMAKE_CURRENT_SOURCE_DIR}/Manual/source/_static/images/SherpaLogo.jpg
  
               DESTINATION  ${CMAKE_INSTALL_DOCDIR}/ COMPONENT doc  )


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
              
