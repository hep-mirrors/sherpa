pkglib_LTLIBRARIES = libAmegic.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

AMEGICFLAGS   =	# -pedantic -Wall -Winline
DEFS          =

libAmegic_la_SOURCES = Helicity.C Topology.C Polarisation.C Pol_Info.C Point.C \
		       Process_Info.C Process_Base.C Single_Process.C \
		       Single_Process_MHV.C Process_Group.C \
		       All_Processes.C Full_Decay_Table.C All_Decays.C Amegic.C
INCLUDES = -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Phasespace \
	   -I../Model \
	   -I../String \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../EXTRA_XS-1.0/Two2Two \
	   -I../../PHASIC++-1.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../BEAM-1.0/Main \
	   -I../../PDF-1.0/Main \
	   -I../../PDF-1.0/LHAPDF \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
pkginclude_HEADERS = Helicity.H Topology.H Polarisation.H Pol_Info.H Point.H \
		 Process_Info.H Process_Base.H Single_Process.H \
		 Single_Process_MHV.H Process_Group.H \
		 All_Processes.H Full_Decay_Table.H All_Decays.H Amegic.H 

EXTRA_DIST    = $(CVSINCLUDE) 
AM_CXXFLAGS   =	$(AMEGICFLAGS)
