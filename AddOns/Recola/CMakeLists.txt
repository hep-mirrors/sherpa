add_library(SherpaRecola SHARED Recola_Interface.C Recola_Born.C  Recola_Virtual.C)
install(TARGETS SherpaRecola DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaRecola       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaRecola SOVERSION 0.0.0)
