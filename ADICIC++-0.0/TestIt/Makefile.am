MAKE    = make
SUBDIRS = data

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
ROOTLDADD     =	-L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTAUXLIBS@
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

if CVS_INCLUDE
CVSINCLUDE = CVS
endif



AM_CXXFLAGS = -pedantic -Wall -Winline -ftemplate-depth-100



bin_PROGRAMS = Flavtest Emissiontest Chaintest Cascadetest

if HDECAY_SUPPORT
HDECAYLIB = -L@libdir@ -lHdecay 
endif
if ISAJET_SUPPORT
ISAJETLIB = -L@libdir@ -lIsajet 
endif
MODELLIBS = -L@libdir@ -lModelMain $(HDECAYLIB) $(ISAJETLIB) 



Flavtest_SOURCES = Flavtest.C
Flavtest_LDADD = @libdir@/libAdicicTool.la \
                 @libdir@/libAdicicInput.la \
                 @libdir@/libAdicicStruct.la \
                 @libdir@/libAdicicHdl.la \
                 @libdir@/libToolsPhys.la \
                 @libdir@/libToolsMath.la \
                 @libdir@/libToolsOrg.la \
                 -lg2c -lm -ldl $(ROOTLDADD) $(MODELLIBS)
#Flavtest_LDFLAGS = SHERPA_LDFLAGS



#-lCint -lCore -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
#-lpthread -lg++
#-L$(rootlibdir)



Emissiontest_SOURCES = Emissiontest.C
Emissiontest_LDADD = @libdir@/libAdicicTool.la \
                     @libdir@/libAdicicInput.la \
                     @libdir@/libAdicicStruct.la \
                     @libdir@/libAdicicHdl.la \
                     @libdir@/libAdicicMain.la \
                     @libdir@/libToolsPhys.la \
                     @libdir@/libToolsMath.la \
                     @libdir@/libToolsOrg.la \
	             -lg2c -lm -ldl $(ROOTLDADD) $(MODELLIBS)



Chaintest_SOURCES = Chaintest.C
Chaintest_LDADD = @libdir@/libAdicicTool.la \
                  @libdir@/libAdicicInput.la \
                  @libdir@/libAdicicStruct.la \
                  @libdir@/libAdicicHdl.la \
                  @libdir@/libAdicicMain.la \
                  @libdir@/libToolsPhys.la \
                  @libdir@/libToolsMath.la \
                  @libdir@/libToolsOrg.la \
	          -lg2c -lm -ldl $(ROOTLDADD) $(MODELLIBS)



Cascadetest_SOURCES = Cascadetest.C
Cascadetest_LDADD = @libdir@/libAdicicTool.la \
                    @libdir@/libAdicicInput.la \
                    @libdir@/libAdicicStruct.la \
                    @libdir@/libAdicicHdl.la \
                    @libdir@/libAdicicMain.la \
                    @libdir@/libToolsPhys.la \
                    @libdir@/libToolsMath.la \
                    @libdir@/libToolsOrg.la \
	            -lg2c -lm -ldl $(ROOTLDADD) $(MODELLIBS)



INCLUDES = -I../Tool \
           -I../Input \
           -I../Structure \
           -I../Handler \
           -I../Main \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org \
           -I../../MODEL-1.0/Main



#include_HEADERS =



EXTRA_DIST = $(CVSINCLUDE)
#EXTRA_DIST = Testrun/Run.dat Testrun/Integration.dat

dist-hook:
	sed -e "s/SHERPA_LDFLAGS/SHERPA_LDFLAGS/1" \
	  < Makefile.am > $(distdir)/Makefile.am


