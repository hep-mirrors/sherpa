rootdir    = $(ROOTSYS)
rootincdir = $(rootdir)/include
rootlibdir = $(rootdir)/lib
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
lib_LTLIBRARIES = libPhasespace.la
libPhasespace_la_SOURCES = Multi_Channel.C Single_Channel.C \
			   ISR_Channel.C FSR_Channel.C Rambo.C Sarge.C \
			   Channel_Basics.C Channel_Elements.C \
			   Phase_Space_Handler.C RamboKK.C \
			   Phase_Space_Generator.C Channel_Generator.C \
			   Phase_Space_Integrator.C 
INCLUDES = -I../ISR \
	   -I../Beam \
	   -I../../AMEGIC++-2.0/Main \
	   -I../../AMEGIC++-2.0/Amplitude \
	   -I../../AMEGIC++-2.0/Amplitude/AmplTools \
	   -I../../AMEGIC++-2.0/Amplitude/Zfunctions \
	   -I../../AMEGIC++-2.0/Model \
	   -I../../AMEGIC++-2.0/String \
	   -I../../EXTRA_XS-1.0/Main \
	   -I../../PDF-1.0/Include \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
           -I../../Analysis-0.7/Main \
           -I../../Analysis-0.7/IDEA \
           -I$(rootincdir) \
           -I$(CLHEPincdir)
noinst_HEADERS = Multi_Channel.H Single_Channel.H \
		 ISR_Channel.H FSR_Channel.H Rambo.H Sarge.H \
		 Channel_Basics.H Channel_Elements.H \
		 Phase_Space_Handler.H RamboKK.H \
		 Phase_Space_Generator.H Channel_Generator.H \
		 Phase_Space_Integrator.H 




