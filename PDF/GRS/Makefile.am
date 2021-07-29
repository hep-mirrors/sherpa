pkglib_LTLIBRARIES    =	libGRSSherpa.la

GITTAG = PDF/GRS
include ../../git.make

AM_FFLAGS=-std=legacy

libGRSSherpa_la_SOURCES       =	\
				Git_Info.C \
				grsg99.f \
				GRSph_Fortran_Interface.C 

localinc_HEADERS      =	\
			GRSph_Fortran_Interface.H 

nobase_dist_pkgdata_DATA         =	\
            GRSGrid/grsg99d.grid \
            GRSGrid/grsg99log.grid \
            GRSGrid/grsg99m.grid \
            GRSGrid/grsg99nlo.grid