
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Org/Makefile.am
 
 #pkglib_LTLIBRARIES = libToolsOrg.la)
 if (SHERPA_ENABLE_GZIP)
 set(GZIPSTREAMSOURCES  Gzip_Stream.C)
 set(GZIPSTREAMHEADERS  Gzip_Stream.H)
 else()
 set(GZIPEXTRADIST  Gzip_Stream.C Gzip_Stream.H)
 endif()
 set(MD5_EXCLUDE  CXXFLAGS_PACKAGES.H CXXFLAGS.H)
 set(GITTAG  ATOOLS/Org)
 #include ../../git.make)
set(Org_SOURCES  MyStrStream.C Return_Value.C Shell_Tools.C STL_Tools.C RUsage.C Getter_Function.C My_MPI.C Exception.C Terminator_Objects.C libzippp.cpp My_File.C File_IO_Base.C Read_Write_Base.C Data_Reader.C
 Yaml_Reader.C Data_Writer.C Run_Parameter.C MyTiming.C Message.C Stacktrace.C IO_Handler.C Info_Key.C Integration_Info.C Library_Loader.C Command_Line_Interface.C Settings.C Scoped_Settings.C Settings_Keys.C
 Settings_Writer.C ${GZIPSTREAMSOURCES} binreloc.c)
set(Org_HEADERS CXXFLAGS_PACKAGES.H CXXFLAGS.H Git_Info.H Return_Value.H Shell_Tools.H STL_Tools.H RUsage.H Getter_Function.H Getter_Function.C My_MPI.H Exception.H Terminator_Objects.H libzippp.h My_File.H
 File_IO_Base.H Read_Write_Base.H Data_Reader.H Yaml_Reader.H Data_Writer.H My_Limits.H Run_Parameter.H Message.H Stacktrace.H MyTiming.H IO_Handler.H MyStrStream.H Info_Key.H Info_Key.inl.H Integration_Info.H
 Integration_Info.inl.H Library_Loader.H Command_Line_Interface.H Command_Line_Options.H Option_Parser.H Settings.H Scoped_Settings.H Settings_Keys.H Settings_Writer.H Strings.H Enum_Flags.H Node.H Node.C
 ${GZIPSTREAMHEADERS} binreloc.h)
 #EXTRA_DIST = $(GZIPEXTRADIST) CXXFLAGS.H.in)
 #dist-hook:)
 #rm -f $(distdir)/CXXFLAGS.H)
 #rm -f $(distdir)/CXXFLAGS_PACKAGES.H)
creategitinfo(Org ../Org/)
  list(TRANSFORM Org_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Org/")
 
set(Org_esources )
if (SHERPA_BUILD_STATIC)
add_library(ToolsOrg_static STATIC ${ToolsOrg_esources}
                             ${Org_SOURCES} 
)
target_include_directories(ToolsOrg_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ToolsOrg_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ToolsOrg_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ToolsOrg_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ToolsOrg)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ToolsOrg SHARED ${ToolsOrg_esources}
                             ${Org_SOURCES} 
)
target_include_directories(ToolsOrg PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ToolsOrg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ToolsOrg DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ToolsOrg        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ToolsOrg SOVERSION 0.0.0)
endif()

if (SHERPA_ENABLE_GZIP)
target_link_libraries(ToolsOrg PRIVATE libzip::zip)
endif()