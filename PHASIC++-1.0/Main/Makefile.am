rootdir    = $(ROOTSYS)
rootincdir = $(rootdir)/include
rootlibdir = $(rootdir)/lib
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
lib_LTLIBRARIES = libPhasespace.la
libPhasespace_la_SOURCES = Phase_Space_Handler.C Phase_Space_Integrator.C \
			   Multi_Channel.C Single_Channel.C \
			   ISR_Channel.C FSR_Channel.C Rambo.C RamboKK.C Sarge.C \
			   Channel_Basics.C Channel_Elements.C 
INCLUDES = -I../ISR \
	   -I../Beam \
	   -I../../AMEGIC++-2.0/Main \
	   -I../../AMEGIC++-2.0/Amplitude \
	   -I../../AMEGIC++-2.0/Amplitude/AmplTools \
	   -I../../AMEGIC++-2.0/Amplitude/Zfunctions \
	   -I../../AMEGIC++-2.0/Phasespace \
	   -I../../AMEGIC++-2.0/Model \
	   -I../../AMEGIC++-2.0/String \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../PDF-1.0/Main \
	   -I../../PDF-1.0/LHAPDF \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
           -I../../Analysis-0.7/Main \
           -I../../Analysis-0.7/IDEA \
           -I$(rootincdir) \
           -I$(CLHEPincdir)
noinst_HEADERS = Phase_Space_Handler.H Phase_Space_Integrator.H \
		 Multi_Channel.H Single_Channel.H \
		 ISR_Channel.H FSR_Channel.H Rambo.H RamboKK.H Sarge.H \
		 Channel_Basics.H Channel_Elements.H 




