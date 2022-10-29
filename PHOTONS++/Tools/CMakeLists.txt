
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Tools/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhotonsTools.la)
 set(GITTAG  PHOTONS++/Tools)
 #include ../../git.make)
set(Tools_SOURCES  Dipole_FF.C Dipole_FI.C Dress_Blob_Base.C Generate_One_Photon.C Weight_Dipole.C Weight_Higher_Order_Corrections.C Weight_Jacobian.C Weight_YFS.C YFS_Form_Factor.C)
set(Tools_HEADERS Dipole_FF.H Dipole_FI.H Dress_Blob_Base.H Generate_One_Photon.H Weight_Base.H Weight_Dipole.H Weight_Higher_Order_Corrections.H Weight_Jacobian.H Weight_YFS.H YFS_Form_Factor.H)
creategitinfo(Tools ../Tools/)
  list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
 
set(Tools_esources )
if (SHERPA_BUILD_STATIC)
add_library(PhotonsTools_static STATIC ${PhotonsTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(PhotonsTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhotonsTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhotonsTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhotonsTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PhotonsTools)
set_target_properties(PhotonsTools_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(PhotonsTools SHARED ${PhotonsTools_esources}
                             ${Tools_SOURCES} 
)
target_include_directories(PhotonsTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhotonsTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhotonsTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhotonsTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhotonsTools SOVERSION 0.0.0)
set_target_properties(PhotonsTools        PROPERTIES DEFINE_SYMBOL "")
endif()
