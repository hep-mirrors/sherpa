##### Automatically or semiautomaticaly generated, do not edit.#########
#Showers/Makefile.am
 
 #pkglib_LTLIBRARIES = libMCatNLOShowers.la)
 set(GITTAG  MCATNLO/Showers)
set(Showers_SOURCES  Splitting_Function_Base.C Splitting_Function_Group.C SF_Coupling.C SF_Lorentz.C Sudakov.C Shower.C Kinematics_Base.C)
set(Showers_HEADERS Splitting_Function_Base.H Splitting_Function_Group.H SF_Coupling.H SF_Lorentz.H SF_Key.H Sudakov.H Shower.H Kinematics_Base.H)
creategitinfo(Showers ../Showers/)
list(TRANSFORM Showers_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Showers/")
list(APPEND Showers_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(MCatNLOShowers SHARED ${Showers_SOURCES} )
target_include_directories(MCatNLOShowers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(MCatNLOShowers)
install(TARGETS MCatNLOShowers DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(MCatNLOShowers        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME MCatNLOShowers SOVERSION 0.0.0)
set_target_properties(MCatNLOShowers        PROPERTIES DEFINE_SYMBOL "")
