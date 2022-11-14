##### Automatically or semiautomaticaly generated, do not edit.#########
#Phasespace/Makefile.am
 
 #pkglib_LTLIBRARIES = libAmegicPSGen.la)
 set(GITTAG  AMEGIC++/Phasespace)
set(Phasespace_SOURCES  Phase_Space_Generator.C Channel_Generator_Base.C Channel_Generator.C Channel_Generator_NPV.C Channel_Generator3V.C Channel_Generator3_NPV.C Channel_Generator_Decays.C Channel_Generator_KK.C
 Channel_Generator_UniV.C)
set(Phasespace_HEADERS Phase_Space_Generator.H Channel_Generator_Base.H Channel_Generator.H Channel_Generator_NPV.H Channel_Generator3V.H Channel_Generator3_NPV.H Channel_Generator_Decays.H Channel_Generator_KK.H
 Channel_Generator_UniV.H)
creategitinfo(Phasespace ../Phasespace/)
list(TRANSFORM Phasespace_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Phasespace/")
list(APPEND Phasespace_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AmegicPSGen SHARED ${Phasespace_SOURCES} )
target_include_directories(AmegicPSGen PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(AmegicPSGen)
install(TARGETS AmegicPSGen DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AmegicPSGen        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AmegicPSGen SOVERSION 0.0.0)
set_target_properties(AmegicPSGen        PROPERTIES DEFINE_SYMBOL "")
