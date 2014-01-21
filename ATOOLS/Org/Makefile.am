pkglib_LTLIBRARIES = libToolsOrg.la

if GZIP_SUPPORT
GZIPSTREAMSOURCES = Gzip_Stream.C
GZIPSTREAMHEADERS = Gzip_Stream.H
else
GZIPEXTRADIST = Gzip_Stream.C Gzip_Stream.H
endif

MD5_EXCLUDE = CXXFLAGS_PACKAGES.H CXXFLAGS.H

SVNTAG = ATOOLS/Org
include ../../svn.make

libToolsOrg_la_SOURCES        = \
				SVN_Info.C \
				MyStrStream.C \
				Return_Value.C \
				Shell_Tools.C \
				STL_Tools.C \
				RUsage.C \
				Getter_Function.C \
				My_MPI.C \
				Exception.C \
				Exception_Handler.C \
				My_File.C \
				File_IO_Base.C \
				Read_Write_Base.C \
				Data_Reader.C \
				Data_Writer.C \
				Run_Parameter.C \
				MyTiming.C \
				Message.C \
				IO_Handler.C \
				Info_Key.C \
				Integration_Info.C \
				Library_Loader.C \
				$(GZIPSTREAMSOURCES) \
				binreloc.c

localinc_HEADERS         =	\
			CXXFLAGS_PACKAGES.H \
			CXXFLAGS.H \
			SVN_Info.H \
			Return_Value.H \
			Smart_Pointer.H \
			Smart_Pointer.C \
			Shell_Tools.H \
			STL_Tools.H \
			RUsage.H \
			Getter_Function.H \
			My_MPI.H \
			Exception.H \
			Exception_Handler.H \
			My_File.H \
			File_IO_Base.H \
			Read_Write_Base.H \
			Read_Write_Base.inl.H \
			Data_Reader.H \
			Data_Writer.H \
			My_Limits.H \
			Run_Parameter.H \
			Message.H \
			MyTiming.H \
			IO_Handler.H \
			MyStrStream.H \
			Info_Key.H \
			Info_Key.inl.H \
			Integration_Info.H \
			Integration_Info.inl.H \
			Library_Loader.H \
			Node.H \
			Node.C \
			$(GZIPSTREAMHEADERS) \
			binreloc.h

libToolsOrg_la_LIBADD = @CONDITIONAL_GZIPLIBS@ @SQLITE3_LDFLAGS@

libToolsOrg_la_CPPFLAGS   =	$(AM_CPPFLAGS) @SQLITE3_CFLAGS@

EXTRA_DIST    = $(GZIPEXTRADIST) CXXFLAGS.H.in

dist-hook:
	rm -f $(distdir)/CXXFLAGS.H
	rm -f $(distdir)/CXXFLAGS_PACKAGES.H
