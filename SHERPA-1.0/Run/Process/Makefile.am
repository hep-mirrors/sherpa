SHERPASYS     =	../../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

SUBDIRS = Dummy
lib_LTLIBRARIES = libProcess.la
libProcess_la_SOURCES =	\
			Set_Channel.C \
			Set_Values.C 
EXTRA_DIST    =	$(CVSINCLUDE) Set_Channel.C.clean Set_Values.C.clean
AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      = $(AMEGICINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) 









