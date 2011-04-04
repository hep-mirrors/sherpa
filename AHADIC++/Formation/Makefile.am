pkglib_LTLIBRARIES = libAhadicFormation.la

SVNTAG = AHADIC++/Formation
include ../../svn.make

libAhadicFormation_la_SOURCES = \
                                SVN_Info.C \
				Cluster_Formation_Handler.C \
                                Gluon_Decayer.C \
                                Cluster_Former.C \
                                Colour_Reconnections.C

localincdir = $(pkgincludedir)/AHADIC++/Formation

localinc_HEADERS   = \
		     Cluster_Formation_Handler.H \
                     Gluon_Decayer.H \
                     Cluster_Former.H \
                     Colour_Reconnections.H

EXTRA_DIST    = @SVNINCLUDE@
