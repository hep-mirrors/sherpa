rootdir    = $(ROOTSYS)
rootincdir = $(rootdir)/include
rootlibdir = $(rootdir)/lib
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
bin_PROGRAMS = Amegic
Amegic_SOURCES = main.C
Amegic_LDADD = @libdir@/libProcess.la \
	       @libdir@/libAmegic.la \
	       @libdir@/libAmplitude.la \
	       @libdir@/libZfunctions.la \
	       @libdir@/libAmplTools.la \
	       @libdir@/libModel.la \
	       @libdir@/libIsajet.la \
	       @libdir@/libString.la \
	       @libdir@/libToolsPhys.la \
	       @libdir@/libToolsMath.la \
	       @libdir@/libToolsOrg.la \
	       @libdir@/libLHAPDF.la \
	       @libdir@/libPDF.la \
	       @libdir@/libExtraXS.la \
	       @libdir@/libExtraXS2_2.la \
	       @libdir@/libPhasespace.la \
	       @libdir@/libISR.la \
	       @libdir@/libBeam.la \
	       -lg2c -lm 
#	       -L$(rootlibdir) -lCint -lCore -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
Amegic_LDFLAGS =  -rdynamic

INCLUDES = -I../Main \
	   -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Model \
	   -I../String \
	   -I../../PHASIC++-1.0/Main \
	   -I../../PHASIC++-1.0/ISR \
	   -I../../PHASIC++-1.0/Beam \
	   -I../../PDF-1.0/Main \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
SUBDIRS = Process
EXTRA_DIST = Testrun/Run.dat Testrun/Integration.dat \
	Testrun/Particle.dat Testrun/Const.dat Testrun/Processes.dat \
	Testrun/ISR.dat Testrun/Hadron.dat Testrun/Isajet.dat \
	Testrun/Run.dat Testrun/Selector.dat Testrun/Selector.dat.sample makeclean 



