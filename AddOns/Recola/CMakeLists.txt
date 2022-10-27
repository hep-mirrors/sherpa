#add_definitions(-DOPENLOOPS_PREFIX="")
#find_package(Recola REQUIRED)

#  message(STATUS "CASCADE: OPENLOOPS_LIBRARIES=${OPENLOOPS_LIBRARIES} OPENLOOPS_PREFIX=${OPENLOOPS_PREFIX} ")



add_library(SherpaRecola SHARED Recola_Interface.C 
				Recola_Born.C 
				Recola_Virtual.C)
#target_link_libraries()
install(TARGETS SherpaRecola DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaRecola       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaRecola SOVERSION 0.0.0)
