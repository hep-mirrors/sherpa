set(Tools_SOURCES Output_Base.C
                  Event_Reader_Base.C
                  Output_LHEF.C
                  Userhook_Base.C
                  Analysis_Interface.C)
if (HEPMC2_REQUIRED)
  list(APPEND Tools_SOURCES HepMC2_Interface.C)
endif()
if (HEPMC3_REQUIRED)
  list(APPEND Tools_SOURCES HepMC3_Interface.C)
endif()

set(Tools_HEADERS Definitions.H
                  Output_Base.H
                  Event_Reader_Base.H
                  HepMC2_Interface.H
                  HepMC3_Interface.H
                  Output_LHEF.H
                  Userhook_Base.H
                  Analysis_Interface.H)

set(foo_SOURCES  Output_Base.C Event_Reader_Base.C HepMC2_Interface.C HepMC3_Interface.C Output_LHEF.C Userhook_Base.C Analysis_Interface.C)
set(foo_HEADERS Definitions.H Output_Base.H Event_Reader_Base.H HepMC2_Interface.H HepMC3_Interface.H Output_LHEF.H Userhook_Base.H Analysis_Interface.H)

sherpa_create_git_info(foo)
 
add_library(SherpaTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaTools SherpaTools)

if (HEPMC2_REQUIRED)
  target_link_libraries(SherpaTools PRIVATE ${HEPMC2_LIBRARIES})
  target_include_directories(SherpaTools PRIVATE ${HEPMC2_INCLUDE_DIRS})
endif()
if (HEPMC3_REQUIRED)
  target_link_libraries(SherpaTools PRIVATE ${HEPMC3_LIB})
  target_include_directories(SherpaTools PRIVATE ${HEPMC3_INCLUDE_DIR})
  if (SHERPA_ENABLE_HEPMC3_ROOT)
    target_link_libraries(SherpaTools PRIVATE ${HEPMC3_ROOTIO_LIB})
  endif()
endif()

