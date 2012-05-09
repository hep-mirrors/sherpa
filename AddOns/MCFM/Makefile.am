
MCFM_SOURCES  =	\
		MCFM_qqb_QQb.C \
		MCFM_gg_h.C \
		MCFM_gg_hgg.C \
		MCFM_qqb_vh.C \
		MCFM_qqb_vv.C \
		MCFM_qqb_v.C \
		MCFM_wbf.C \
		MCFM_gb_hb.C \
		MCFM_qqb_wgam_v.C \
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

EXTRA_DIST = @SVNINCLUDE@ $(MCFM_EXTRA_DIST) install_mcfm.sh
