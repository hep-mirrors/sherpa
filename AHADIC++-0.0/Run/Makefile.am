AHADICSYS  = ../..

AHADICFLAGS   =	-pedantic -Wall -Winline 

ATOOLSDIR     =	$(AHADICSYS)/ATOOLS-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0

AHADICINCS    =	-I../Main -I../ClusterFormation -I../ClusterDecays
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PHASICINCS    =	-I$(AHADICSYS)/PHASIC++-1.0/Main
PDFINCS       =	-I$(AHADICSYS)/PDF-1.0/Main
BEAMINCS      =	-I$(AHADICSYS)/BEAM-1.0/Main
MODELINCS     =	-I$(AHADICSYS)/MODEL-1.0/Main
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers

SYSLIBS	      =	-lg2c -lm -ldl
ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 
PHASICLIBS    =	-L@libdir@ -lPhasespace 
AHADICLIBS    =	-L@libdir@ -lAhadicMain -lClusterFormation -lClusterDecays
PDFLIBS       =	-L@libdir@ -lDUPDF -lSudakov -lPDF -lLHAPDF \
                           -lCTEQ -lMRST -lGRV -lRemnant
BEAMLIBS      =	-L@libdir@ -lBeam 
MODELLIBS     =	-L@libdir@ -lModelMain -lIsajet -lHdecay 
APACICLIBS    =	-L@libdir@ -lApacicShowers -lApacicMain 

bin_PROGRAMS  =	Cluster
Cluster_SOURCES = Main.C 
Cluster_LDADD  = $(AHADICLIBS) $(PHASICLIBS) $(ATOOLSLIBS) $(SYSLIBS) \
		 $(MODELLIBS) $(PDFLIBS) $(BEAMLIBS) $(APACICLIBS)
Cluster_LDFLAGS        =	-rdynamic
AM_CXXFLAGS   =	$(AHADICFLAGS)
INCLUDES      = $(AHADICINCS) $(PHASICINCS) $(ATOOLSINCS) $(MODELINCS) \
		$(PDFINCS) $(BEAMINCS) $(APACICINCS)
