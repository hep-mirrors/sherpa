pkglib_LTLIBRARIES = libAhadicMain.la

SVNTAG = AHADIC++__Main
include ../../svn.make

libAhadicMain_la_SOURCES = \
                           SVN_Info.C \
			   Ahadic.C

localincdir = $(pkgincludedir)/AHADIC++/Main

localinc_HEADERS   = \
                     SVN_Info.C.in \
		     Ahadic.H

EXTRA_DIST    = @SVNINCLUDE@
