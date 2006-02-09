SHERPASYS  = ../..
CSSFLAGS   =	-pedantic -Wall -Winline 

CSSDIR        =	$(SHERPASYS)/CS_SHOWER++-0.0
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0

ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     = -I$(SHERPASYS)/MODEL-1.0/Main
CSSINCS       = -I$(CSSDIR)/Tools -I$(CSSDIR)/Showers -I$(CSSDIR)/Main

SYSLIBS	      =	-lg2c -lm -ldl

ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 
CSSLIBS       =	-L@libdir@ -lCSSTools -lCSShowers -lCSMain
MODELLIBS     =	-L@libdir@ -lModelMain -lIsajet -lHdecay 

bin_PROGRAMS    = CSSTest
CSSTest_SOURCES = Main.C 
CSSTest_LDADD   = $(CSSLIBS) $(ATOOLSLIBS) $(MODELLIBS) $(SYSLIBS)
CSSTest_LDFLAGS = -rdynamic
AM_CXXFLAGS   =	$(CSSFLAGS) 
INCLUDES      = $(CSSINCS) $(ATOOLSINCS) $(MODELINCS) 
EXTRA_DIST    =	Particle.dat 

dist-hook:
	sed -e "s/-rdynamic/SHERPA_LDFLAGS/1" \
	  < Makefile.am > $(distdir)/Makefile.am
