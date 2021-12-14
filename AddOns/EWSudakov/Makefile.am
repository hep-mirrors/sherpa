if USING__EWSud

SUBDIRS = Tests

endif

GITTAG = AddOns/EWSudakov
include ../../git.make

EWSUD_SOURCES = \
        Variation_Generator.C \
        Variation_Weight_Names.C \
	Coefficient_Checker.C \
	KFactor_Checker.C \
	Comix_Interface.C \
	HE_Comix_Interface.C \
	EWGroupConstants.C \
	EWSudakov.C \
	EWSudakov_Amplitudes.C \
	Git_Info.C \
	KFactor.C \
	EWSudakov_Calculator.C \
	Clustered_EWSudakov_Calculator.C
EWSUD_ADDS = \
        Variation_Generator.H \
        Variation_Weight_Names.H \
	Coefficient_Checker.H \
	KFactor_Checker.H \
	Comix_Interface.H \
	HE_Comix_Interface.H \
	EWGroupConstants.H \
	EWSudakov.H \
	EWSudakov_Amplitudes.H \
	KFactor.H \
	EWSudakov_Calculator.H \
	Clustered_EWSudakov_Calculator.H

if USING__EWSud

pkglib_LTLIBRARIES    =	libSherpaEWSud.la
libSherpaEWSud_la_SOURCES  =	$(EWSUD_SOURCES)
localinc_HEADERS      = $(EWSUD_ADDS)
EWSUD_EXTRA_DIST = $(EWSUD_ADDS)

endif

EXTRA_DIST = $(EWSUD_EXTRA_DIST)
