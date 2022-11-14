##### Automatically or semiautomaticaly generated, do not edit.#########
#Cross_Sections/Makefile.am
 
 set(GITTAG  SHRiMPS/Cross_Sections)
set(Cross_Sections_SOURCES  Cross_Sections.C Sigma_Base.C Sigma_Total.C Sigma_Inelastic.C Sigma_Elastic.C Sigma_SD.C Sigma_Partonic.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Cross_Sections)
set(Cross_Sections_HEADERS Cross_Sections.H Sigma_Base.H Sigma_Total.H Sigma_Inelastic.H Sigma_Elastic.H Sigma_SD.H Sigma_Partonic.H)
sherpa_create_git_info(Cross_Sections ../Cross_Sections/)
list(TRANSFORM Cross_Sections_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Cross_Sections/")
list(APPEND Cross_Sections_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsXsecs SHARED ${Cross_Sections_SOURCES} )
sherpa_handle_shared_library(ShrimpsXsecs ShrimpsXsecs)
