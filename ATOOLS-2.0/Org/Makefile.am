SHERPASYS = ../..

ATOOLSFLAGS   =	-DUSING__COLOUR # -pedantic -Wall -Winline
DEFS          =

ATOOLSINCS    =	-I../Math -I../Phys
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
ROOTFLAGS     = -Wno-long-long
endif

lib_LTLIBRARIES = libToolsOrg.la
libToolsOrg_la_SOURCES        = \
				Shell_Tools.C \
				Object.C \
				Getter_Function.C \
				Exception.C \
				Exception_Handler.C \
				My_File.C \
				File_IO_Base.C \
				Read_Write_Base.C \
				Data_Reader.C \
				Data_Writer.C \
				Data_Return.C \
				Data_Read.C \
				Data_Collector.C \
				Run_Parameter.C \
				MyTiming.C \
				Message.C \
		 		Sort.C \
				prof.cc \
				IO_Handler.C \
				Info_Key.C \
				Integration_Info.C \
				Debugger.C \
				My_Root.C
pkginclude_HEADERS       =	\
			CXXFLAGS.H \
			Smart_Pointer.H \
			Shell_Tools.H \
			Object.H \
			Getter_Function.H \
			Exception.H \
			Exception_Handler.H \
			My_File.H \
			File_IO_Base.H \
			Read_Write_Base.H \
			Read_Write_Base.inl.H \
			Data_Reader.H \
			Data_Writer.H \
			Data_Read.H \
			My_Limits.H \
			Data_Return.H \
			Data_Collector.H \
			Run_Parameter.H \
			Message.H \
			List_Algorithms.H \
			Sort.H \
			MyTiming.H \
			prof.hh \
			prof_exampl.cc \
			IO_Handler.H \
			MyStrStream.H \
			Info_Key.H \
			Info_Key.inl.H \
			Integration_Info.H \
			Integration_Info.inl.H \
			Debugger.H \
			Debugger.inl.H \
			Node.H \
			Node.C \
			My_Root.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ATOOLSFLAGS) $(ROOTFLAGS) \
		-DSHERPA_INCLUDE_PATH=\"$(pkgincludedir)\" \
		-DSHERPA_BINARY_PATH=\"$(bindir)\" \
		-DSHERPA_PDFS_PATH=\"$(pkgdatadir)\"
INCLUDES      = $(MODELINCS) $(ATOOLSINCS) $(ROOTINCS)
