##### Automatically or semiautomaticaly generated, do not edit.#########
#Electron/Makefile.am
 
 #pkglib_LTLIBRARIES = libPDFESherpa.la)
 set(GITTAG  PDF/Electron)
set(Electron_SOURCES  PDF_Electron.C)
set(Electron_HEADERS PDF_Electron.H)
sherpa_create_git_info(Electron ../Electron/)
list(TRANSFORM Electron_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Electron/")
list(APPEND Electron_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PDFESherpa SHARED ${Electron_SOURCES} )
sherpa_handle_shared_library(PDFESherpa PDFESherpa)
