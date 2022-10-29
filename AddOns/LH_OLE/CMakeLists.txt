add_library(SherpaLHOLE SHARED LH_OLE_Interface.C 
                        LH_OLE_Communicator.C)
install(TARGETS SherpaLHOLE DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaLHOLE       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaLHOLE SOVERSION 0.0.0)
set_target_properties(SherpaLHOLE       PROPERTIES DEFINE_SYMBOL "")
