#Manually edited
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#HZTool/Makefile.am
 
 set(GITTAG  AddOns/HZTool)
 #include ../../git.make)
set(HZTool_SOURCES  hznorm.F hztool_sherpa.F HZTool_Interface.C)
set(HZTOOL_ADDS   HZTool_Wrapper.H)
 #AM_FFLAGS = -fno-automatic @CONDITIONAL_HZTOOLINCS@)
 if (SHERPA_ENABLE_HZTOOL)
 if (SHERPA_ENABLE_CERNLIB)
 #pkglib_LTLIBRARIES = libSherpaHZToolAnalysis.la)
set(HZTool_SOURCES ${HZTool_SOURCES})
 #HZTOOL_#EXTRA_DIST = $(HZTOOL_ADDS)
 else()
 #HZTOOL_#EXTRA_DIST = $(HZTOOL_ADDS) ${HZTOOL_SOURCES})
 endif()
 else()
 #HZTOOL_#EXTRA_DIST = $(HZTOOL_ADDS) ${HZTOOL_SOURCES})
 endif()
 #EXTRA_DIST = $(#HZTOOL_#EXTRA_DIST)
creategitinfo(HZTool ../HZTool/)
list(TRANSFORM HZTool_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../HZTool/")
list(APPEND HZTool_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaHZToolAnalysis SHARED 
                             ${HZTool_SOURCES} 
)
target_include_directories(SherpaHZToolAnalysis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaHZToolAnalysis DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaHZToolAnalysis        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaHZToolAnalysis SOVERSION 0.0.0)
set_target_properties(SherpaHZToolAnalysis        PROPERTIES DEFINE_SYMBOL "")
target_include_directories( SherpaHZToolAnalysis PRIVATE ${HZTOOL_INCLUDE_DIR} )
target_link_libraries( SherpaHZToolAnalysis PRIVATE ${HZTOOL_LIBRARY} 
${CERN_kernlib_DYNAMIC_LIBRARY}
${CERN_mathlib_DYNAMIC_LIBRARY}
${CERN_packlib_DYNAMIC_LIBRARY}
${MOTIF_LIBRARIES}
)
