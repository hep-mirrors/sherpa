SHERPASYS = ../..

ATOOLSFLAGS   =	-DUSING__COLOUR # -pedantic -Wall -Winline
DEFS          =

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
ROOTINCS      =	-I@ROOTINCDIR@
ROOTLDADD     =	-L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTAUXLIBS@
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

bin_PROGRAMS  =	PI_Demo Calc Color
PI_Demo_SOURCES       =	PI_Demo.C
PI_Demo_LDADD         =	$(ATOOLSLIBS) $(SYSLIBS) $(ROOTLDADD) 
PI_Demo_LDFLAGS       =	@ldflags@
Color_SOURCES =	Color.C
Color_LDADD   =	$(ATOOLSLIBS) $(SYSLIBS) $(ROOTLDADD) 
Color_LDFLAGS =	@ldflags@
Calc_SOURCES  =	Calc.C
Calc_LDADD    =	$(ATOOLSLIBS) $(SYSLIBS) $(ROOTLDADD)
Calc_LDFLAGS  =	@ldflags@
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(MODELFLAGS) $(ATOOLSFLAGS) $(ROOTFLAGS)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS) $(ROOTINCS)

