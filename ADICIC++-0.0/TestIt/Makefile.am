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



bin_PROGRAMS = Flavtest



Flavtest_SOURCES = Flavtest.C
Flavtest_LDADD = @libdir@/libAdicicTool.la \
                 @libdir@/libAdicicStruct.la \
                 @libdir@/libToolsPhys.la \
                 @libdir@/libToolsMath.la \
                 @libdir@/libToolsOrg.la \
                 -lg2c -lm
#Flavtest_LDFLAGS = -rdynamic



#-lCint -lCore -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -ldl
#-lpthread -lg++
#-L$(rootlibdir)



INCLUDES = -I../Tool \
           -I../Structure \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org



#noinst_HEADERS =



EXTRA_DIST = $(CVSINCLUDE)
#EXTRA_DIST = Testrun/Run.dat Testrun/Integration.dat







