find_package(OpenLoops REQUIRED)
message(STATUS "SHERPA: OPENLOOPS_LIBRARIES=${OPENLOOPS_LIBRARIES} OPENLOOPS_PREFIX=${OPENLOOPS_PREFIX} ")
add_library(SherpaOpenLoops SHARED OpenLoops_Virtual.C 
            OpenLoops_Born.C 
            OpenLoops_Interface.C 
            Color_Correlated_ME2.C 
            Spin_Color_Correlated_ME2.C 
            GGH_Process_Manager.C 
            GGH_KFactor_Setter.C 
            GGH_Decay_KFactor_Setter.C)

install(TARGETS SherpaOpenLoops DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaOpenLoops       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaOpenLoops SOVERSION 0.0.0)

target_link_libraries(SherpaOpenLoops  PRIVATE ${OPENLOOPS_LIBRARY})
