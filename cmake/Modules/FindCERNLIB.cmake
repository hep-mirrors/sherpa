 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 #     Modified by Andrii Verbytskyi, MPI fuer Physik, 2016
 ################################################################################
# - Try to find CERNLIB
# Once done this will define
#
#  CERNLIB_FOUND - system has CERNLIB
#### (not needed)  CERNLIB_INCLUDE_DIR - the CERNLIB include directory
#  CERNLIB_LIBRARIES - The libraries needed to use CERNLIB
#### (not needed)  CERNLIB_DEFINITIONS - Compiler switches required for using CERNLIB
#

if (CERNLIB_INCLUDE_DIR AND CERNLIB_LIBRARY_DIR)
  SET (CERNLIB_INCLUDE_DIR CERNLIB_INCLUDE_DIR-NOTFOUND)
  SET (CERNLIB_LIB_DIR CERNLIB_LIB_DIR-NOTFOUND)
  SET (CERNLIB_PLISTS_LIB_DIR CERNLIB_PLISTS_LIB_DIR-NOTFOUND)
endif (CERNLIB_INCLUDE_DIR AND CERNLIB_LIBRARY_DIR)



FIND_PATH(CERNLIB_INCLUDE_DIR NAMES CERNLIB PATHS  
  $ENV{CERN_ROOT}/include  
  /usr/include/cernlib/2006
  /usr/include/cernlib/2022
  NO_DEFAULT_PATH 
)

FIND_PATH(CERNLIB_LIBRARY_DIR NAMES  libpacklib.a   libkernlib.a libmathlib.a     libpacklib.so   libkernlib.so libmathlib.so PATHS  
  $ENV{CERN_ROOT}/lib
  /usr/lib64/cernlib/2006/lib
  /usr/lib64/cernlib/2022/lib
  /usr/lib/cernlib/2006/lib
  /usr/lib/cernlib/2022/lib
  NO_DEFAULT_PATH
)

if (CERNLIB_LIBRARY_DIR)
   set(CERNLIB_FOUND TRUE)
endif (CERNLIB_LIBRARY_DIR)

set(cernlibs kernlib packlib mathlib graflib grafX11)  
if (CERNLIB_FOUND)


set(CERN_DYNAMIC_LIBRARIES)
set(CERN_STATIC_LIBRARIES)
foreach(_cpt ${cernlibs})

  find_library(CERN_${_cpt}_DYNAMIC_LIBRARY lib${_cpt}.so ${_cpt} HINTS ${CERNLIB_LIBRARY_DIR})
  if(CERN_${_cpt}_DYNAMIC_LIBRARY)
    mark_as_advanced(CERN_${_cpt}_DYNAMIC_LIBRARY)
    list(APPEND CERN_DYNAMIC_LIBRARIES ${CERN_${_cpt}_DYNAMIC_LIBRARY})
  endif()

  find_library(CERN_${_cpt}_STATIC_LIBRARY lib${_cpt}.a ${_cpt} HINTS ${CERNLIB_LIBRARY_DIR})
  if(CERN_${_cpt}_STATIC_LIBRARY)
    mark_as_advanced(CERN_${_cpt}_STATIC_LIBRARY)
    list(APPEND CERN_STATIC_LIBRARIES ${CERN_${_cpt}_STATIC_LIBRARY})
  endif()
 
endforeach()

list(REMOVE_DUPLICATES CERN_DYNAMIC_LIBRARIES)
list(REMOVE_DUPLICATES CERN_STATIC_LIBRARIES)
  
else (CERNLIB_FOUND)
  if (CERNLIB_FIND_REQUIRED)
  message(FATAL_ERROR "Looking for CERNLIB... - Not found")
  endif (CERNLIB_FIND_REQUIRED)
endif (CERNLIB_FOUND)


SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${CERNLIB_LIBRARY_DIR})


