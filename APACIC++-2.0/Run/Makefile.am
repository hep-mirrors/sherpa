bin_PROGRAMS = Apacic
rootdir    = $(ROOTSYS)
rootincdir = $(rootdir)/include
rootlibdir = $(rootdir)/lib
CLHEPdir   = /usr/local/CLHEP
#  was ../../AEXTERN before !!!
CLHEPlibdir   = $(CLHEPdir)/lib
CLHEPincdir   = $(CLHEPdir)/include
Apacic_SOURCES = Apacic.C Main.C
Apacic_LDADD =  \
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
	 -lg2c -lm \
	 -L$(rootlibdir) -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
Apacic_LDFLAGS =  -rdynamic
INCLUDES = -I../Main \
	   -I../Hard_Processes \
	   -I../Hard_Processes/Process_Library \
	   -I../Hard_Processes/Integration_Tools \
	   -I../Hard_Interface \
	   -I../Showers \
	   -I../Run \
	   -I../Soft_Interface \
	   -I../../PDF-1.0/Include \
	   -I../../EXTRA_XS-1.0/Main \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org \
	   -I../../AMEGIC++-2.0/Amplitude \
	   -I../../AMEGIC++-2.0/Amplitude/AmplTools \
	   -I../../AMEGIC++-2.0/Amplitude/Zfunctions \
	   -I../../AMEGIC++-2.0/Main \
	   -I../../AMEGIC++-2.0/Model \
	   -I../../AMEGIC++-2.0/String \
	   -I../../PHASIC++-1.0/Main \
	   -I../../PHASIC++-1.0/ISR \
	   -I../../PHASIC++-1.0/Beam \
	   -I../../Analysis/Main \
	   -I../../Analysis/IDEA \
	   -I$(rootincdir) -I$(CLHEPincdir)
noinst_HEADERS = Apacic.H
EXTRA_DIST =  Particle.dat Hadron.dat Const.dat Run.dat Integration.dat ISR.dat 
# cor01.dat
#	     cor01.dat cor02.dat cor03.dat cor04.dat cor05.dat cor06.dat 
#	     cor07.dat cor08.dat cor09.dat cor10.dat cor11.dat cor12.dat






