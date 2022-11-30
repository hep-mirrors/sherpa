##### Automatically or semiautomaticaly generated, do not edit.#########
#CT14/Makefile.am
 
set(GITTAG  PDF/CT14)

set(CT14_SOURCES CT14Pdf.f
                 CT14_Fortran_Interface.C)

sherpa_create_git_info(CT14)
 
add_library(CT14Sherpa SHARED ${CT14_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CT14Sherpa CT14Sherpa)
