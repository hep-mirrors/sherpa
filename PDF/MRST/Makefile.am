pkglib_LTLIBRARIES    = \
			libMRST99Sherpa.la \
			libMRST01LOSherpa.la \
			libMRST04QEDSherpa.la

SVNTAG = PDF__MRST99
include ../../svn.make

libMRST99Sherpa_la_SOURCES    =	\
				SVN_Info.C \
				mrst99.cc \
				PDF_MRST99.C
libMRST01LOSherpa_la_SOURCES  =	\
				SVN_Info.C \
				mrst2001lo.f \
				PDF_MRST01LO.C
libMRST04QEDSherpa_la_SOURCES =	\
				SVN_Info.C \
				mrst2004qed.f \
				PDF_MRST04QED.C

localinc_HEADERS      = \
			SVN_Info.C.in \
			mrst99.h \
			PDF_MRST99.H \
			PDF_MRST01LO.H \
			PDF_MRST04QED.H

nobase_dist_pkgdata_DATA         =	\
			MRST99Grid/cor01.dat \
			MRST99Grid/cor02.dat \
			MRST99Grid/cor03.dat \
			MRST99Grid/cor04.dat \
			MRST99Grid/cor05.dat \
			MRST99Grid/cor06.dat \
			MRST99Grid/cor07.dat \
			MRST99Grid/cor08.dat \
			MRST99Grid/cor09.dat \
			MRST99Grid/cor10.dat \
			MRST99Grid/cor11.dat \
			MRST99Grid/cor12.dat \
			MRST01Grid/lo2002.dat \
			MRST04Grid/qed6-10gridn.dat \
			MRST04Grid/qed6-10gridp.dat
