pkglib_LTLIBRARIES       =	libLundTools.la

PYTHIASOURCES = pythia-6.4.18.f \
		pyplst.f \
		Lund_Interface.C \
		Lund_Decay_Handler.C
PYTHIAHEADERS = Lund_Wrapper.H Lund_Interface.H Lund_Decay_Handler.H

MD5_EXCLUDE   =	pythia-6.4.18.f

GITTAG = SHERPA/LundTools
include ../../git.make

if PYTHIA_SUPPORT
libLundTools_la_SOURCES       =	$(PYTHIASOURCES) Git_Info.C
localinc_HEADERS        =      $(PYTHIAHEADERS) Dummy.C
else
libLundTools_la_SOURCES       = Git_Info.C Dummy.C
localinc_HEADERS        =       $(PYTHIASOURCES) $(PYTHIAHEADERS)
endif
