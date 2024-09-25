set(Tools_SOURCES Output_Base.C
                  Event_Reader_Base.C
                  Output_LHEF.C
                  Userhook_Base.C
                  Analysis_Interface.C)
if (HepMC3_FOUND)
  list(APPEND Tools_SOURCES HepMC3_Interface.C)
endif()

set(Tools_HEADERS Definitions.H
                  Output_Base.H
                  Event_Reader_Base.H
                  HepMC3_Interface.H
                  Output_LHEF.H
                  Userhook_Base.H
                  Analysis_Interface.H)

set(foo_SOURCES  Output_Base.C Event_Reader_Base.C HepMC3_Interface.C Output_LHEF.C Userhook_Base.C Analysis_Interface.C)
set(foo_HEADERS Definitions.H Output_Base.H Event_Reader_Base.H HepMC3_Interface.H Output_LHEF.H Userhook_Base.H Analysis_Interface.H)

sherpa_create_git_info(foo)
 
add_library(SherpaTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaTools SherpaTools)

if (HepMC3_FOUND)
  target_link_libraries(SherpaTools PRIVATE HepMC3::HepMC3 ${HepMC3_rootIO})
endif()

