rootdir    = $(ROOTSYS)
rootincdir = $(rootdir)/include
rootlibdir = $(rootdir)/lib
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
bin_PROGRAMS   = Mocaic
Mocaic_SOURCES = main.C
Mocaic_LDADD =  \
         @libdir@/libMocaic.la \
         @libdir@/libMEPSInterface.la \
         @libdir@/libHadInterface.la \
	 @libdir@/libMain.la \
	 @libdir@/libAmegic.la \
	 @libdir@/libAmplitude.la \
	 @libdir@/libPhasespace.la \
	 @libdir@/libZfunctions.la \
	 @libdir@/libAmplTools.la \
	 @libdir@/libModel.la \
	 @libdir@/libIsajet.la \
	 @libdir@/libProcess.la \
	 @libdir@/libString.la \
	 @libdir@/libISR.la \
	 @libdir@/libBeam.la \
	 @libdir@/libHardProc.la \
	 @libdir@/libHardInterface.la \
	 @libdir@/libShowers.la \
	 @libdir@/libPDF.la \
	 @libdir@/libExtraXS.la \
	 @libdir@/libExtraXS2_2.la \
	 @libdir@/libToolsPhys.la \
	 @libdir@/libToolsMath.la \
	 @libdir@/libToolsOrg.la \
	 -lg2c -lm 
#	 -L$(rootlibdir) -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
Mocaic_LDFLAGS =  -rdynamic
INCLUDES = -I../Mocaic/ \
	   -I../ME_PS_Interface/ \
	   -I../Had_Interface/ \
	   -I../../AMEGIC++-2.0/Amplitude \
	   -I../../AMEGIC++-2.0/Amplitude/AmplTools \
	   -I../../AMEGIC++-2.0/Amplitude/Zfunctions \
	   -I../../AMEGIC++-2.0/Main \
	   -I../../AMEGIC++-2.0/Model \
	   -I../../AMEGIC++-2.0/String \
	   -I../../APACIC++-2.0/Main \
	   -I../../APACIC++-2.0/Hard_Interface \
	   -I../../APACIC++-2.0/Hard_Processes \
	   -I../../APACIC++-2.0/Showers \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../PDF-1.0/Include \
	   -I../../PHASIC++-1.0/Main \
	   -I../../PHASIC++-1.0/ISR \
	   -I../../PHASIC++-1.0/Beam \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org 
SUBDIRS    = Process
EXTRA_DIST = Integration.dat Processes.dat Selector.dat \
	     Particle.dat Hadron.dat Const.dat Run.dat ISR.dat cor01.dat makeclean

