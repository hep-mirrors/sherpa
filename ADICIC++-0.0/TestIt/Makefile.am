MAKE    = make
SUBDIRS = data

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTLIBS      =	-L$(ROOTSYS)/lib -lCint -lCore -lHist -lGraf \
		-lGraf3d -lGpad -lTree -lRint -lPostscript \
		-lMatrix -lPhysics -lGui -ldl 
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

if CVS_INCLUDE
CVSINCLUDE = CVS
endif



AM_CXXFLAGS = -pedantic -Wall -Winline



bin_PROGRAMS = Flavtest Emissiontest Chaintest Cascadetest



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
                 -lg2c -lm -ldl $(ROOTLIBS)
#Flavtest_LDFLAGS = -rdynamic



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
                     @libdir@/libModelMain.la \
                     @libdir@/libHdecay.la \
                     @libdir@/libIsajet.la \
                     -lg2c -lm -ldl $(ROOTLIBS)



Chaintest_SOURCES = Chaintest.C
Chaintest_LDADD = @libdir@/libAdicicTool.la \
                  @libdir@/libAdicicInput.la \
                  @libdir@/libAdicicStruct.la \
                  @libdir@/libAdicicHdl.la \
                  @libdir@/libAdicicMain.la \
                  @libdir@/libToolsPhys.la \
                  @libdir@/libToolsMath.la \
                  @libdir@/libToolsOrg.la \
                  @libdir@/libModelMain.la \
                  @libdir@/libHdecay.la \
                  @libdir@/libIsajet.la \
                  -lg2c -lm -ldl $(ROOTLIBS)



Cascadetest_SOURCES = Cascadetest.C
Cascadetest_LDADD = @libdir@/libAdicicTool.la \
                    @libdir@/libAdicicInput.la \
                    @libdir@/libAdicicStruct.la \
                    @libdir@/libAdicicHdl.la \
                    @libdir@/libAdicicMain.la \
                    @libdir@/libToolsPhys.la \
                    @libdir@/libToolsMath.la \
                    @libdir@/libToolsOrg.la \
                    @libdir@/libModelMain.la \
                    @libdir@/libHdecay.la \
                    @libdir@/libIsajet.la \
                    -lg2c -lm -ldl $(ROOTLIBS)



INCLUDES = -I../Tool \
           -I../Input \
           -I../Structure \
           -I../Handler \
           -I../Main \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org \
           -I../../MODEL-1.0/Main



#noinst_HEADERS =



EXTRA_DIST = $(CVSINCLUDE)
#EXTRA_DIST = Testrun/Run.dat Testrun/Integration.dat






