
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#NNLO/Makefile.am
 
 set(GITTAG  AddOns/NNLO)
 #include ../../git.make)
 #pkglib_LTLIBRARIES = libNNLOqT.la)
set(NNLO_SOURCES  polylog.c Tools.C QT_Selector.C DIS_Selector.C Scale_Setter.C coeffqt.C param.C hard.C DY_KFactor.C higgsfullsm.C H_KFactor.C DISinclu.C DIS_KFactor.C DIS_Scale.C)
#NNLOqT_#EXTRA_DIST = Tools.H QT_Selector.H coeffqt.H param.H hard.H higgsfullsm.H DISinclu.H
 #EXTRA_DIST = $(#NNLOqT_#EXTRA_DIST)
creategitinfo(NNLO ../NNLO/)
list(TRANSFORM NNLO_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../NNLO/")
list(APPEND NNLO_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(NNLOqT SHARED 
                             ${NNLO_SOURCES} 
)
target_include_directories(NNLOqT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(NNLOqT)
install(TARGETS NNLOqT DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(NNLOqT        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME NNLOqT SOVERSION 0.0.0)
set_target_properties(NNLOqT        PROPERTIES DEFINE_SYMBOL "")
