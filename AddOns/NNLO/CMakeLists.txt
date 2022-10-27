#add_definitions(-DOPENLOOPS_PREFIX="")
#find_package(HZToolAnalysis REQUIRED)

#  message(STATUS "CASCADE: OPENLOOPS_LIBRARIES=${OPENLOOPS_LIBRARIES} OPENLOOPS_PREFIX=${OPENLOOPS_PREFIX} ")



add_library(NNLOqT SHARED polylog.c 
	Tools.C 
	QT_Selector.C 
	DIS_Selector.C 
	Scale_Setter.C 
	coeffqt.C 
	param.C 
	hard.C 
        DY_KFactor.C 
	higgsfullsm.C 
	H_KFactor.C 
	DISinclu.C 
	DIS_KFactor.C 
	DIS_Scale.C)

install(TARGETS NNLOqT DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(NNLOqT       PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME NNLOqT SOVERSION 0.0.0)
