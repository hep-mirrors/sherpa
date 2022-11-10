
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Perturbative/Makefile.am
 
 #pkglib_LTLIBRARIES = libAmisicPerturbative.la)
 set(GITTAG  AMISIC++/Perturbative)
 #include ../../git.make)
set(Perturbative_SOURCES  Single_Collision_Handler.C MI_Processes.C MI_Process_Group.C MI_Process.C QCD_Processes.C QED_Processes.C)
set(Perturbative_HEADERS Single_Collision_Handler.H MI_Processes.H MI_Process_Group.H MI_Process.H QCD_Processes.H QED_Processes.H)
creategitinfo(Perturbative ../Perturbative/)
list(TRANSFORM Perturbative_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Perturbative/")
list(APPEND Perturbative_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AmisicPerturbative SHARED 
                             ${Perturbative_SOURCES} 
)
target_include_directories(AmisicPerturbative PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(AmisicPerturbative)
install(TARGETS AmisicPerturbative DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AmisicPerturbative        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AmisicPerturbative SOVERSION 0.0.0)
set_target_properties(AmisicPerturbative        PROPERTIES DEFINE_SYMBOL "")
