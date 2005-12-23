lib_LTLIBRARIES = libAmegicPSGen.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

AMEGICFLAGS   =	# -pedantic -Wall #-Winline
DEFS          =

libAmegicPSGen_la_SOURCES = Phase_Space_Generator.C Channel_Generator_Base.C \
			    Channel_Generator.C Channel_Generator_NPV.C \
                            Channel_Generator3V.C Channel_Generator3_NPV.C \
	 		    Channel_Generator_Decays.C Channel_Generator_KK.C

INCLUDES = -I../Main \
	   -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Model \
	   -I../String \
	   -I../../PHASIC++-1.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../BEAM-1.0/Main \
	   -I../../PDF-1.0/Main \
	   -I../../PDF-1.0/LHAPDF \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../EXTRA_XS-1.0/Two2Two \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
include_HEADERS = Phase_Space_Generator.H Channel_Generator_Base.H \
		 Channel_Generator.H Channel_Generator_NPV.H \
		 Channel_Generator3V.H Channel_Generator3_NPV.H \
		 Channel_Generator_Decays.H Channel_Generator_KK.H

EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMEGICFLAGS)
