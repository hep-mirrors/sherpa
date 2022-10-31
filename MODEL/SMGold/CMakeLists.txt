
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#SMGold/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaSMGold.la)
 set(GITTAG  MODEL/SMGold)
 #include ../../git.make)
set(SMGold_SOURCES  Model.C Identity_2_1.C T_3_2_1.C VSS1.C VSS2.C VVVV2.C VVVV3.C VVVV4.C)
creategitinfo(SMGold ../SMGold/)
list(TRANSFORM SMGold_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SMGold/")
list(APPEND SMGold_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(SMGold_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaSMGold SHARED ${SherpaSMGold_esources}
                             ${SMGold_SOURCES} 
)
target_include_directories(SherpaSMGold PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSMGold PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSMGold DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSMGold        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSMGold SOVERSION 0.0.0)
set_target_properties(SherpaSMGold        PROPERTIES DEFINE_SYMBOL "")
endif()
