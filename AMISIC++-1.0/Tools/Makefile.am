SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0

lib_LTLIBRARIES               =	libAmisicTools.la
libAmisicTools_la_SOURCES     = \
				Read_Write_Base.C \
				Data_Reader.C \
				Data_Writer.C \
				Variable.C 
CXXFLAGS      =	-pedantic -Wall -Winline
INCLUDES      = \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org
noinst_HEADERS        =	\
			My_String.H \
			My_Functional.H \
			My_Map.H \
			My_Vector.H \
			My_IO_Stream.H \
			My_Stringstream.H \
			My_Filestream.H \
			Type.H \
			Type.IC \
			Read_Write_Base.H \
			Read_Write_Base.IC \
			Data_Reader.H \
			Data_Reader.IC \
			Data_Writer.H \
			Variable.H \
			Variable.IC \
			Scaling.H \
			Scaling.C \
			Axis.H \
			Axis.IC \
			Axis.C \
			Data_To_Function.H \
			Data_To_Function.IC \
			Data_To_Function.C 
