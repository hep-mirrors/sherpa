add_library(SherpaGoSam SHARED 
        GoSam_Virtual.C 
                        GoSam_Interface.C)

target_include_directories(SherpaGoSam PRIVATE ${PROJECT_BINARY_DIR})
install(TARGETS SherpaGoSam DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaGoSam      PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaGoSam SOVERSION 0.0.0)
set_target_properties(SherpaGoSam       PROPERTIES DEFINE_SYMBOL "")
