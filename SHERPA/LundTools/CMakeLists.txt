#Manually edited
##### Automatically or semiautomaticaly generated, do not edit.#########

set(PYTHIASOURCES   pyplst.f Lund_Interface.C)
 set(PYTHIAHEADERS  Lund_Wrapper.H Lund_Interface.H)
 set(GITTAG  SHERPA/LundTools)
 if (SHERPA_ENABLE_PYTHIA6)
  set(LundTools_SOURCES ${PYTHIASOURCES} )
  set(LundTools_HEADERS ${PYTHIAHEADERS} Dummy.C)
  sherpa_create_git_info(LundTools ../LundTools/)
  list(TRANSFORM LundTools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LundTools/")
  list(APPEND LundTools_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/../LundTools/pythia-6.4.18.f")
 else()
  set(LundTools_SOURCES  Dummy.C)
  set(LundTools_HEADERS ${PYTHIASOURCES} ${PYTHIAHEADERS})
  sherpa_create_git_info(LundTools ../LundTools/)
  list(TRANSFORM LundTools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LundTools/")
 endif()
list(APPEND LundTools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(LundTools SHARED ${LundTools_SOURCES} )
sherpa_handle_shared_library(LundTools LundTools)

