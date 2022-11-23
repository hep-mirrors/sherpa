##### Automatically or semiautomaticaly generated, do not edit.#########
#Main/Makefile.am
 
set(GITTAG  AMEGIC++/Main)
set(Main_SOURCES  ColorSC.C Helicity.C Topology.C Polarisation.C Pol_Info.C Point.C Process_Tags.C Amegic_Base.C Process_Base.C Single_Process.C Single_Process_MHV.C Single_Process_External.C Single_Process_Combined.C
 Process_Group.C Amegic.C)
set(Main_HEADERS Tools.H ColorSC.H Helicity.H Topology.H Polarisation.H Pol_Info.H Point.H Process_Tags.H Amegic_Base.H Process_Base.H Single_Process.H Single_Process_MHV.H Single_Process_External.H Single_Process_Combined.H
 Process_Group.H)
sherpa_create_git_info(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(Amegic SHARED ${Main_SOURCES} )
sherpa_handle_shared_library(Amegic Amegic)
