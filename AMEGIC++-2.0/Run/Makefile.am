rootdir    = $(ROOTSYS)
rootincdir = $(rootdir)/include
rootlibdir = $(rootdir)/lib
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
bin_PROGRAMS = Amegic

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

Amegic_SOURCES = main.C Environment.C
Amegic_LDADD = @libdir@/libProcess.la \
	       @libdir@/libAmegic.la \
	       @libdir@/libAmegicPSGen.la \
	       @libdir@/libAmplitude.la \
	       @libdir@/libZfunctions.la \
	       @libdir@/libAmplTools.la \
	       @libdir@/libModel.la \
	       @libdir@/libIsajet.la \
	       @libdir@/libString.la \
	       @libdir@/libToolsPhys.la \
	       @libdir@/libToolsMath.la \
	       @libdir@/libToolsOrg.la \
	       @libdir@/libModelMain.la \
	       @libdir@/libHdecay.la \
	       @libdir@/libIsajet.la \
	       @libdir@/libModel.la \
	       @libdir@/libDUPDF.la \
	       @libdir@/libApacicMain.la \
	       @libdir@/libApacicShowers.la \
	       @libdir@/libPDF.la \
	       @libdir@/libSudakov.la \
	       @libdir@/libLHAPDF.la \
	       @libdir@/libPhasespace.la \
	       @libdir@/libBeam.la \
	       -lg2c -lm -ldl
#	       -L$(rootlibdir) -lCint -lCore -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
Amegic_LDFLAGS =  -rdynamic
noinst_HEADERS = Environment.H 

INCLUDES = -I../Main \
	   -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Phasespace \
	   -I../Model \
	   -I../String \
	   -I../../PHASIC++-1.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../MODEL-1.0/Hdecay \
	   -I../../MODEL-1.0/Isajet \
	   -I../../BEAM-1.0/Main \
	   -I../../PDF-1.0/Main \
	   -I../../PDF-1.0/LHAPDF \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
SUBDIRS = Process
EXTRA_DIST = makeclean $(CVSINCLUDE)



