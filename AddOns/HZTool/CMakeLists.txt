#add_definitions(-DOPENLOOPS_PREFIX="")
#find_package(HZToolAnalysis REQUIRED)

#  message(STATUS "CASCADE: OPENLOOPS_LIBRARIES=${OPENLOOPS_LIBRARIES} OPENLOOPS_PREFIX=${OPENLOOPS_PREFIX} ")



add_library(SherpaHZToolAnalysis SHARED hznorm.F 
			hztool_sherpa.F 
			HZTool_Interface.C)
target_include_directories( SherpaHZToolAnalysis PRIVATE /usr/include/hztool )
install(TARGETS SherpaHZToolAnalysis DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaHZToolAnalysis       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaHZToolAnalysis SOVERSION 0.0.0)
