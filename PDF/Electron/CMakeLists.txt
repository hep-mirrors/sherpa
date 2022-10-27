
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Electron/Makefile.am
 
 #pkglib_LTLIBRARIES = libPDFESherpa.la)
 set(GITTAG  PDF/Electron)
 #include ../../git.make)
set(Electron_SOURCES  PDF_Electron.C)
set(Electron_HEADERS PDF_Electron.H)
creategitinfo(Electron ../Electron/)
  list(TRANSFORM Electron_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Electron/")
 
set(Electron_esources )
if (SHERPA_BUILD_STATIC)
add_library(PDFESherpa_static STATIC ${PDFESherpa_esources}
                             ${Electron_SOURCES} 
)
target_include_directories(PDFESherpa_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PDFESherpa_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PDFESherpa_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PDFESherpa_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PDFESherpa)
endif()
if (SHERPA_BUILD_SHARED)
add_library(PDFESherpa SHARED ${PDFESherpa_esources}
                             ${Electron_SOURCES} 
)
target_include_directories(PDFESherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PDFESherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PDFESherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PDFESherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PDFESherpa SOVERSION 0.0.0)
endif()
