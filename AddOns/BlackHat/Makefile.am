SVNTAG = AddOns/BlackHat
include ../../svn.make

BLACKHAT_SOURCES      =	\
			SVN_Info.C \
			BlackHat_Tree.C \
			BlackHat_Virtual.C \
			BlackHat_Interface.C

BLACKHAT_ADDS         =	\
			BlackHat_Tree.H \
			BlackHat_Virtual.H

if BLACKHAT_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaBlackHat.la libSherpaWhiteHat.la

libSherpaBlackHat_la_SOURCES  =	$(BLACKHAT_SOURCES)
libSherpaWhiteHat_la_SOURCES  =	$(BLACKHAT_SOURCES)

libSherpaBlackHat_la_LIBADD   =	@CONDITIONAL_BLACKHATLIBS@
libSherpaWhiteHat_la_LIBADD   =	@CONDITIONAL_BLACKHATLIBS@

libSherpaBlackHat_la_CPPFLAGS =	$(AM_CPPFLAGS) @CONDITIONAL_BLACKHATINCS@
libSherpaWhiteHat_la_CPPFLAGS =	$(AM_CPPFLAGS) @CONDITIONAL_BLACKHATINCS@

BLACKHAT_EXTRA_DIST = $(BLACKHAT_ADDS)
else
BLACKHAT_EXTRA_DIST = $(BLACKHAT_SOURCES) $(BLACKHAT_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(BLACKHAT_EXTRA_DIST)
