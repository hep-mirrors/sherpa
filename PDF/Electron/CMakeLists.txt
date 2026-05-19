set(Electron_SOURCES PDF_Electron.C)

set(Electron_HEADERS PDF_Electron.H)

sherpa_create_git_info(Electron ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(PDFESherpa SHARED ${Electron_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(PDFESherpa PDFESherpa)
add_dependencies(PDFESherpa Electron_Git_Info)
