rootdir    = $(ROOTSYS)
rootincdir = $(rootdir)/include
rootlibdir = $(rootdir)/lib
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
bin_PROGRAMS   = Sherpa
Sherpa_SOURCES = Main.C
Sherpa_LDADD =  \
         @libdir@/libSherpaMain.la \
         @libdir@/libSherpaTools.la \
         @libdir@/libSherpaInitialization.la \
         @libdir@/libSherpaPerturbativePhysics.la \
         @libdir@/libSherpaSoftPhysics.la \
         @libdir@/libLundInterface.la \
         @libdir@/libSherpaSingleEvents.la \
	 @libdir@/libModelMain.la \
	 @libdir@/libHdecay.la \
	 @libdir@/libIsajet.la \
	 @libdir@/libAmegic.la \
	 @libdir@/libAmplitude.la \
	 @libdir@/libAmegicPSGen.la \
	 @libdir@/libPhasespace.la \
	 @libdir@/libZfunctions.la \
	 @libdir@/libAmplTools.la \
	 @libdir@/libModel.la \
	 @libdir@/libProcess.la \
	 @libdir@/libString.la \
	 @libdir@/libApacicShowers.la \
	 @libdir@/libApacicMain.la \
	 @libdir@/libBeam.la \
	 @libdir@/libLHAPDF.la \
	 @libdir@/libPDF.la \
	 @libdir@/libExtraXS.la \
	 @libdir@/libExtraXS2_2.la \
	 @libdir@/libToolsPhys.la \
	 @libdir@/libToolsMath.la \
	 @libdir@/libToolsOrg.la \
	 -lg2c -lm 
#	 -L$(rootlibdir) -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
Sherpa_LDFLAGS =  -rdynamic
INCLUDES = -I../../SHERPA-1.0/Main \
	   -I../../SHERPA-1.0/Tools \
	   -I../../SHERPA-1.0/PerturbativePhysics \
           -I../../SHERPA-1.0/SoftPhysics \
           -I../../SHERPA-1.0/SoftPhysics/Lund \
	   -I../../SHERPA-1.0/Initialization \
	   -I../../SHERPA-1.0/Single_Events \
	   -I../../AMEGIC++-2.0/Amplitude \
	   -I../../AMEGIC++-2.0/Amplitude/AmplTools \
	   -I../../AMEGIC++-2.0/Amplitude/Zfunctions \
           -I../../AMEGIC++-2.0/Phasespace \
	   -I../../AMEGIC++-2.0/Main \
	   -I../../AMEGIC++-2.0/Model \
	   -I../../AMEGIC++-2.0/String \
	   -I../../APACIC++-2.0/Main \
	   -I../../APACIC++-2.0/Hard_Interface \
	   -I../../APACIC++-2.0/Showers \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../BEAM-1.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../MODEL-1.0/Hdecay \
	   -I../../PDF-1.0/Main \
	   -I../../PDF-1.0/LHAPDF \
	   -I../../PHASIC++-1.0/Main \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org 
SUBDIRS    = Process
EXTRA_DIST = Integration.dat Processes.dat Selector.dat Model.dat Beam.dat \
	     Particle.dat Hadron.dat Const.dat Run.dat ISR.dat ME.dat makeclean

