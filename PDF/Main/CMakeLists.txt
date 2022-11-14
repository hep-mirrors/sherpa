##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libPDF.la)
 set(GITTAG  PDF/Main)
set(Main_SOURCES  PDF_Base.C ISR_Base.C Intact.C Structure_Function.C ISR_Handler.C Jet_Criterion.C Cluster_Definitions_Base.C NLOMC_Base.C Shower_Base.C Photon_PDF_Base.C)
set(Main_HEADERS PDF_Base.H ISR_Base.H Intact.H Structure_Function.H ISR_Handler.H Jet_Criterion.H Cluster_Definitions_Base.H NLOMC_Base.H Shower_Base.H Photon_PDF_Base.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PDF SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(PDF PDF)
