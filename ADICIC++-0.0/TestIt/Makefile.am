MAKE    = make
SUBDIRS = data

rootdir     = $(ROOTSYS)
rootincdir  = $(rootdir)/include
rootlibdir  = $(rootdir)/lib
CLHEPlibdir = $(CLHEPDIR)/lib
CLHEPincdir = $(CLHEPDIR)/include



if CVS_INCLUDE
CVSINCLUDE = CVS
endif



AM_CXXFLAGS = -pedantic -Wall -Winline



bin_PROGRAMS = Flavtest Emissiontest Chaintest



Flavtest_SOURCES = Flavtest.C
Flavtest_LDADD = @libdir@/libAdicicTool.la \
                 @libdir@/libAdicicInput.la \
                 @libdir@/libAdicicStruct.la \
                 @libdir@/libAdicicHdl.la \
                 @libdir@/libToolsPhys.la \
                 @libdir@/libToolsMath.la \
                 @libdir@/libToolsOrg.la \
                 @libdir@/libModelMain.la \
                 @libdir@/libHdecay.la \
                 @libdir@/libIsajet.la \
                 -lg2c -lm -ldl
#Flavtest_LDFLAGS = -rdynamic



#-lCint -lCore -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
#-lpthread -lg++
#-L$(rootlibdir)



Emissiontest_SOURCES = Emissiontest.C
Emissiontest_LDADD = @libdir@/libAdicicTool.la \
                     @libdir@/libAdicicInput.la \
                     @libdir@/libAdicicStruct.la \
                     @libdir@/libAdicicHdl.la \
                     @libdir@/libToolsPhys.la \
                     @libdir@/libToolsMath.la \
                     @libdir@/libToolsOrg.la \
                     @libdir@/libModelMain.la \
                     @libdir@/libHdecay.la \
                     @libdir@/libIsajet.la \
                     -lg2c -lm -ldl



Chaintest_SOURCES = Chaintest.C
Chaintest_LDADD = @libdir@/libAdicicTool.la \
                  @libdir@/libAdicicInput.la \
                  @libdir@/libAdicicStruct.la \
                  @libdir@/libAdicicHdl.la \
                  @libdir@/libToolsPhys.la \
                  @libdir@/libToolsMath.la \
                  @libdir@/libToolsOrg.la \
                  @libdir@/libModelMain.la \
                  @libdir@/libHdecay.la \
                  @libdir@/libIsajet.la \
                  -lg2c -lm -ldl



INCLUDES = -I../Tool \
           -I../Input \
           -I../Structure \
           -I../Handler \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org \
           -I../../MODEL-1.0/Main



#noinst_HEADERS =



EXTRA_DIST = $(CVSINCLUDE)
#EXTRA_DIST = Testrun/Run.dat Testrun/Integration.dat





