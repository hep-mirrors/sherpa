add_library(SherpaHZToolAnalysis SHARED hznorm.F 
            hztool_sherpa.F 
            HZTool_Interface.C)
target_include_directories( SherpaHZToolAnalysis PRIVATE /usr/include/hztool ) #!FIXME
install(TARGETS SherpaHZToolAnalysis DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaHZToolAnalysis       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaHZToolAnalysis SOVERSION 0.0.0)
