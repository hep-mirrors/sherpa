find_package(BlackHat REQUIRED)
message(STATUS "SHERPA: BLACKHAT_LIBRARIES=${BLACKHAT_LIBRARIES} BLACKHAT_INCLUDE_DIRS=${BLACKHAT_INCLUDE_DIRS}")
add_library(SherpaBlackHat SHARED BlackHat_Tree.C 
            BlackHat_Virtual.C 
            BlackHat_Interface.C)
target_link_libraries(SherpaBlackHat PRIVATE ${BLACKHAT_LIBRARIES})
target_include_directories(SherpaBlackHat PRIVATE ${BLACKHAT_INCLUDE_DIRS})
install(TARGETS SherpaBlackHat DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaBlackHat       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaBlackHat SOVERSION 0.0.0)
set_target_properties(SherpaBlackHat       PROPERTIES DEFINE_SYMBOL "")
