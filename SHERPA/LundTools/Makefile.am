pkglib_LTLIBRARIES       =	libLundTools.la

PYTHIASOURCES = pythia-6.4.18.f \
		pyplst.f \
		Lund_Interface.C
PYTHIAHEADERS = Lund_Wrapper.H Lund_Interface.H

SVNTAG = SHERPA__LundTools
include ../../svn.make

if PYTHIA_SUPPORT
libLundTools_la_SOURCES       =	$(PYTHIASOURCES) SVN_Info.C
localinc_HEADERS        =	$(PYTHIAHEADERS) SVN_Info.C.in Dummy.C
else
libLundTools_la_SOURCES       = SVN_Info.C Dummy.C
localinc_HEADERS        =       $(PYTHIASOURCES) $(PYTHIAHEADERS) SVN_Info.C.in
endif
