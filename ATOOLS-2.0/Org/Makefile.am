SHERPASYS = ../..

ATOOLSFLAGS   =	-DUSING__COLOUR -pedantic -Wall -Winline

ATOOLSINCS    =	-I../Math -I../Phys
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

lib_LTLIBRARIES = libToolsOrg.la
libToolsOrg_la_SOURCES        = \
				Shell_Tools.C \
				Object.C \
				Getter_Function.C \
				Exception.C \
				Exception_Handler.C \
				File_IO_Base.C \
				Read_Write_Base.C \
				Data_Reader.C \
				Data_Writer.C \
				Data_Read.C \
				Run_Parameter.C \
				FastList.C \
				MyTiming.C \
				Message.C \
		 		Sort.C \
				prof.cc \
				IO_Handler.C \
				Type.C \
				Info_Key.C \
				Integration_Info.C \
				Debugger.C \
				My_Root.C
noinst_HEADERS        =	\
			Shell_Tools.H \
			Object.H \
			Getter_Function.H \
			Exception.H \
			Exception_Handler.H \
			Exception_Handler.inl.H \
			File_IO_Base.H \
			Read_Write_Base.H \
			Read_Write_Base.inl.H \
			Data_Reader.H \
			Data_Reader.inl.H \
			Data_Writer.H \
			Data_Read.H \
			My_Limits.H \
			Data_Return.H \
			Run_Parameter.H \
			Message.H \
			FastList.H \
			List_Algorithms.H \
			Sort.H \
			MyTiming.H \
			prof.hh \
			prof_exampl.cc \
			IO_Handler.H \
			MyStrStream.H \
			Type.H \
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
AM_CXXFLAGS   =	$(ATOOLSFLAGS) $(ROOTFLAGS)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS) $(ROOTINCS)
