pkglib_LTLIBRARIES = libToolsOrg.la

if GZIP_SUPPORT
GZIPSTREAMSOURCES = Gzip_Stream.C
GZIPSTREAMHEADERS = Gzip_Stream.H
else
GZIPEXTRADIST = Gzip_Stream.C Gzip_Stream.H
endif

MD5_EXCLUDE = CXXFLAGS_PACKAGES.H CXXFLAGS.H

GITTAG = ATOOLS/Org
include ../../git.make

libToolsOrg_la_SOURCES        = \
				Git_Info.C \
				MyStrStream.C \
				Return_Value.C \
				Shell_Tools.C \
				STL_Tools.C \
				RUsage.C \
				Getter_Function.C \
				My_MPI.C \
				Exception.C \
				Terminator_Objects.C \
				libzippp.cpp \
				My_File.C \
				File_IO_Base.C \
				Read_Write_Base.C \
				Data_Reader.C \
				Yaml_Reader.C \
				Data_Writer.C \
				Run_Parameter.C \
				MyTiming.C \
				Message.C \
				Stacktrace.C \
				IO_Handler.C \
				Info_Key.C \
				Integration_Info.C \
				Library_Loader.C \
				Command_Line_Interface.C \
				Settings.C \
				Scoped_Settings.C \
				Settings_Keys.C \
				Settings_Writer.C \
				$(GZIPSTREAMSOURCES) \
				binreloc.c

localinc_HEADERS         =	\
			CXXFLAGS_PACKAGES.H \
			CXXFLAGS.H \
			Git_Info.H \
			Return_Value.H \
			Shell_Tools.H \
			STL_Tools.H \
			RUsage.H \
			Getter_Function.H \
			Getter_Function.C \
			My_MPI.H \
			Exception.H \
			Terminator_Objects.H \
			libzippp.h \
			My_File.H \
			File_IO_Base.H \
			Read_Write_Base.H \
			Data_Reader.H \
			Yaml_Reader.H \
			Data_Writer.H \
			My_Limits.H \
			Run_Parameter.H \
			Message.H \
			Stacktrace.H \
			MyTiming.H \
			IO_Handler.H \
			MyStrStream.H \
			Info_Key.H \
			Info_Key.inl.H \
			Integration_Info.H \
			Integration_Info.inl.H \
			Library_Loader.H \
			Command_Line_Interface.H \
			Command_Line_Options.H \
			Option_Parser.H \
			Settings.H \
			Scoped_Settings.H \
			Settings_Keys.H \
			Settings_Writer.H \
			Enum_Flags.H \
			Node.H \
			Node.C \
			$(GZIPSTREAMHEADERS) \
			binreloc.h

libToolsOrg_la_LIBADD = @CONDITIONAL_GZIPLIBS@ @LIBZIP_LDFLAGS@

libToolsOrg_la_CPPFLAGS   =	$(AM_CPPFLAGS) @CONDITIONAL_GZIPINCS@ @LIBZIP_CPPFLAGS@

EXTRA_DIST    = $(GZIPEXTRADIST) CXXFLAGS.H.in

dist-hook:
	rm -f $(distdir)/CXXFLAGS.H
	rm -f $(distdir)/CXXFLAGS_PACKAGES.H
