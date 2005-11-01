SHERPASYS = ../..

ATOOLSFLAGS   =	-DUSING__COLOUR -pedantic -Wall -Winline

ATOOLSINCS    =	-I../Org -I../Math -I../Phys
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if MODEL_SUPPORT
MODELFLAGS    =
else
MODELFLAGS    = -DUSING__ATOOLS_only
endif

SYSLIBS	      =	-lm -ldl
ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTLIBS      =	-L$(ROOTSYS)/lib -lCint -lCore -lHist -lGraf \
		-lGraf3d -lGpad -lTree -lRint -lPostscript \
		-lMatrix -lPhysics -lGui -ldl 
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

bin_PROGRAMS  =	PI_Demo Calc Color
PI_Demo_SOURCES       =	PI_Demo.C
PI_Demo_LDADD         =	$(ATOOLSLIBS) $(SYSLIBS) $(ROOTLIBS) 
PI_Demo_LDFLAGS       =	-rdynamic
Color_SOURCES =	Color.C
Color_LDADD   =	$(ATOOLSLIBS) $(SYSLIBS) $(ROOTLIBS) 
Color_LDFLAGS =	-rdynamic
Calc_SOURCES  =	Calc.C
Calc_LDADD    =	$(ATOOLSLIBS) $(SYSLIBS) $(ROOTLIBS)
Calc_LDFLAGS  =	-rdynamic
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(MODELFLAGS) $(ATOOLSFLAGS) $(ROOTFLAGS)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS) $(ROOTINCS)

dist-hook:
	sed -e "s/-rdynamic/SHERPA_LDFLAGS/1" \
	  < Makefile.am > $(distdir)/Makefile.am
