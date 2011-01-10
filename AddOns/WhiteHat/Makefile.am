SVNTAG = AddOns__WhiteHat
include ../../svn.make

WHITEHAT_SOURCES      =	\
			SVN_Info.C \
			WhiteHat_Tree.C \
			WhiteHat_Virtual.C \
			WhiteHat_Interface.C

WHITEHAT_ADDS         =	\
			SVN_Info.C.in \
			WhiteHat_Tree.H \
			WhiteHat_Virtual.H

if WHITEHAT_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaWhiteHat.la

libSherpaWhiteHat_la_SOURCES  =	$(WHITEHAT_SOURCES)

libSherpaWhiteHat_la_LIBADD   =	@CONDITIONAL_WHITEHATLIBS@

libSherpaWhiteHat_la_CPPFLAGS =	$(AM_CPPFLAGS) @CONDITIONAL_WHITEHATINCS@

WHITEHAT_EXTRA_DIST = $(WHITEHAT_ADDS)
else
WHITEHAT_EXTRA_DIST = $(WHITEHAT_SOURCES) $(WHITEHAT_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(WHITEHAT_EXTRA_DIST)
