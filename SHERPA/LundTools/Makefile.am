pkglib_LTLIBRARIES       =	libLundTools.la

PYTHIASOURCES = pythia6214.f pydata.f pyalps.f pypdfu.f pyplst.f pymult.f \
		Lund_Interface.C
PYTHIAHEADERS = Lund_Wrapper.H Lund_Interface.H

if PYTHIA_SUPPORT
libLundTools_la_SOURCES       =	$(PYTHIASOURCES)
localinc_HEADERS        =	$(PYTHIAHEADERS) Dummy.C
else
libLundTools_la_SOURCES       = Dummy.C
localinc_HEADERS        =       $(PYTHIASOURCES) $(PYTHIAHEADERS)
endif
