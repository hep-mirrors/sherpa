pkglib_LTLIBRARIES = libAmegic.la

SVNTAG = AMEGIC++__Main
include ../../svn.make

libAmegic_la_SOURCES = \
			SVN_Info.C \
			AMEGIC_SCT.C ColorSC.C \
			Helicity.C Topology.C Polarisation.C Pol_Info.C Point.C \
			Process_Tags.C Amegic_Base.C Process_Base.C \
			Single_Process.C Single_Process_MHV.C  Single_Process_External.C \
			Process_Group.C Amegic.C

localinc_HEADERS   = \
			AMEGIC_SCT.H ColorSC.H \
			Helicity.H Topology.H Polarisation.H Pol_Info.H Point.H \
			Process_Tags.H Amegic_Base.H Process_Base.H \
			Single_Process.H Single_Process_MHV.H  Single_Process_External.H \
			Process_Group.H 

dist_pkgdata_SCRIPTS  =	makelibs
