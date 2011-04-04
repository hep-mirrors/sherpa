pkglib_LTLIBRARIES = libAhadicDecays.la

SVNTAG = AHADIC++/Decays
include ../../svn.make

libAhadicDecays_la_SOURCES = \
                             SVN_Info.C \
			     Cluster_Decay_Handler.C \
                             Cluster_Decay_Analysis.C \
                             Cluster_Part.C 

localincdir = $(pkgincludedir)/AHADIC++/Decays

localinc_HEADERS   = \
		     Cluster_Decay_Handler.H \
                     Cluster_Decay_Analysis.H \
                     Cluster_Part.H 

EXTRA_DIST    = @SVNINCLUDE@
