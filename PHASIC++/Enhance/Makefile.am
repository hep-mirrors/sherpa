
pkglib_LTLIBRARIES    =	libPhasicEnhance.la

SVNTAG = PHASIC++/Enhance
include ../../svn.make

libPhasicEnhance_la_SOURCES   =	\
				SVN_Info.C \
				Enhance_Observable_Base.C \
				Variable_Enhance_Observable.C \
				Fastjet_Enhance_Observable.C

localinc_HEADERS      =	\
			Enhance_Observable_Base.H

libPhasicEnhance_la_LIBADD    =	@CONDITIONAL_FASTJETLIBS@

libPhasicEnhance_la_CPPFLAGS  =	$(AM_CPPFLAGS) @CONDITIONAL_FASTJETINCS@
