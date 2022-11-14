##### Automatically or semiautomaticaly generated, do not edit.#########
#Electron/Makefile.am
 
 #pkglib_LTLIBRARIES = libPDFESherpa.la)
 set(GITTAG  PDF/Electron)
set(Electron_SOURCES  PDF_Electron.C)
set(Electron_HEADERS PDF_Electron.H)
creategitinfo(Electron ../Electron/)
list(TRANSFORM Electron_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Electron/")
list(APPEND Electron_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PDFESherpa SHARED ${Electron_SOURCES} )
target_include_directories(PDFESherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(PDFESherpa)
install(TARGETS PDFESherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PDFESherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PDFESherpa SOVERSION 0.0.0)
set_target_properties(PDFESherpa        PROPERTIES DEFINE_SYMBOL "")
