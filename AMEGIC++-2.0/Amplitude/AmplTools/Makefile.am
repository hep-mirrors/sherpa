lib_LTLIBRARIES = libAmplTools.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

AMEGICFLAGS  =	-pedantic -Wall -Winline

libAmplTools_la_SOURCES = CFColor.C
INCLUDES = -I../../Amplitude \
	   -I../../Amplitude/AmplTools \
	   -I../../Amplitude/Zfunctions \
	   -I../../Main \
	   -I../../Model \
	   -I../../String \
	   -I../../ISR \
	   -I../../../MODEL-1.0/Main \
	   -I../../../ATOOLS-2.0/Phys \
	   -I../../../ATOOLS-2.0/Math \
	   -I../../../ATOOLS-2.0/Org 
noinst_HEADERS = CFColor.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMEGICFLAGS)
