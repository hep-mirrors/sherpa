lib_LTLIBRARIES = libAmegicPSGen.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

libAmegicPSGen_la_SOURCES = Phase_Space_Generator.C Channel_Generator.C
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
noinst_HEADERS = Phase_Space_Generator.H Channel_Generator.H

EXTRA_DIST    = $(CVSINCLUDE)