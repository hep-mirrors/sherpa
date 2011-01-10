pkglib_LTLIBRARIES = libToolsOrg.la

if GZIP_SUPPORT
GZIPSTREAMSOURCES = Gzip_Stream.C
GZIPSTREAMHEADERS = Gzip_Stream.H
else
GZIPEXTRADIST = Gzip_Stream.C Gzip_Stream.H
endif

Run_Parameter.C: SVN_Info.H.in
	@if ! which svn > /dev/null; then exit 0; fi; \
	url=$$(svn info | awk '{ if ($$1=="URL:") { \
	  split($$2,a,"svn/sherpa/"); sub("/ATOOLS","",a[2]); \
	  sub("/Org","",a[2]); gsub("/","\\/",a[2]); print a[2]; } }'); \
	sed -e's/SVNBRANCH/'$$url'/1' -e's/SVNREVISION/'$$(svnversion)'/1' \
	  < $< > $(basename $<).tmp; \
	if ! diff $(basename $<).tmp $(basename $<) > /dev/null 2>&1; \
	  then mv $(basename $<).tmp $(basename $<); \
	  else rm $(basename $<).tmp; fi;

.PHONY: SVN_Info.H.in

data_DATA     =	SVN_Info.H

libToolsOrg_la_SOURCES        = \
				MyStrStream.C \
				Return_Value.C \
				Shell_Tools.C \
				STL_Tools.C \
				Getter_Function.C \
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
				binreloc.c \
				My_Root.C

localinc_HEADERS         =	\
			SVN_Info.H \
			CXXFLAGS_PACKAGES.H \
			CXXFLAGS.H \
			Return_Value.H \
			Smart_Pointer.H \
			Smart_Pointer.C \
			Shell_Tools.H \
			STL_Tools.H \
			Getter_Function.H \
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
			List_Algorithms.H \
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
			binreloc.h \
			My_Root.H

libToolsOrg_la_LIBADD = @CONDITIONAL_GZIPLIBS@

libToolsOrg_la_CXXFLAGS = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libToolsOrg_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@ $(GZIPEXTRADIST) CXXFLAGS.H.in

dist-hook:
	rm -f $(distdir)/CXXFLAGS.H
	rm -f $(distdir)/CXXFLAGS_PACKAGES.H
