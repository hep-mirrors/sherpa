find_package(ROOT REQUIRED)
 
 message(STATUS "TMDLIB: ROOT_VERSION=${ROOT_VERSION} ROOT_LIBRARIES=${ROOT_LIBRARIES} ROOT_INCLUDE_DIRS=${ROOT_INCLUDE_DIRS}")

 
add_library(SherpaRootOutput SHARED Output_RootNtuple.C)
target_link_libraries(SherpaRootOutput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaRootOutput PRIVATE ${ROOT_INCLUDE_DIRS})
install(TARGETS SherpaRootOutput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaRootOutput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaRootOutput SOVERSION 0.0.0)



add_library(SherpaERootOutput SHARED Output_RootNtuple.C)
target_link_libraries(SherpaERootOutput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaERootOutput PRIVATE ${ROOT_INCLUDE_DIRS})
install(TARGETS SherpaERootOutput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaERootOutput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaERootOutput SOVERSION 0.0.0)


add_library(SherpaEDRootOutput SHARED Output_RootNtuple.C)
target_link_libraries(SherpaEDRootOutput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaEDRootOutput PRIVATE ${ROOT_INCLUDE_DIRS})
install(TARGETS SherpaEDRootOutput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaEDRootOutput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaEDRootOutput SOVERSION 0.0.0)


add_library(SherpaRootInput SHARED RootNtuple_Reader.C)
target_link_libraries(SherpaRootInput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaRootInput PRIVATE ${ROOT_INCLUDE_DIRS})
install(TARGETS SherpaRootInput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaRootInput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaRootInput SOVERSION 0.0.0)



add_library(SherpaERootInput SHARED RootNtuple_Reader.C)
target_link_libraries(SherpaERootInput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaERootInput PRIVATE ${ROOT_INCLUDE_DIRS})
install(TARGETS SherpaERootInput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaERootInput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaERootInput SOVERSION 0.0.0)


add_library(SherpaEDRootInput SHARED RootNtuple_Reader.C)
target_link_libraries(SherpaEDRootInput PRIVATE ${ROOT_LIBRARIES})
target_include_directories(SherpaEDRootInput PRIVATE ${ROOT_INCLUDE_DIRS})
install(TARGETS SherpaEDRootInput DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaEDRootInput       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaEDRootInput SOVERSION 0.0.0)
