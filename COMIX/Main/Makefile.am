pkglib_LTLIBRARIES = libComix.la

SVNTAG = COMIX/Main
include ../../svn.make

libComix_la_SOURCES   =	\
			SVN_Info.C \
			Process_Base.C \
			Process_Group.C \
			Single_Process.C \
			Single_Dipole_Term.C \
			Comix.C

localinc_HEADERS      =	\
			Process_Base.H \
			Process_Group.H \
			Single_Process.H \
			Single_Dipole_Term.H

dist_pkgdata_SCRIPTS  =	init_nlo.sh
