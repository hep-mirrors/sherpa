SVNTAG = AddOns__MCFM
include ../../svn.make

MCFM_SOURCES  =	\
		SVN_Info.C \
		MCFM_qqb_z_v.C \
		MCFM_qqb_w_v.C \
		MCFM_qqb_w1jet_v.C \
		MCFM_qqb_zz_v.C \
		MCFM_qqb_ww_v.C \
		MCFM_Interface.C

MCFM_ADDS     =	\
		MCFM_Wrapper.H

if MCFM_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaMCFM.la

libSherpaMCFM_la_SOURCES      =	$(MCFM_SOURCES)

libSherpaMCFM_la_LIBADD       =	@CONDITIONAL_MCFMLIBS@

libSherpaMCFM_la_CPPFLAGS     =	$(AM_CPPFLAGS)

MCFM_EXTRA_DIST = $(MCFM_ADDS)
else
MCFM_EXTRA_DIST = $(MCFM_SOURCES) $(MCFM_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(MCFM_EXTRA_DIST)
