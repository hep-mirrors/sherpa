
pkglib_LTLIBRARIES    = libPhasicSelectors.la

SVNTAG = PHASIC++/Selectors
include ../../svn.make

libPhasicSelectors_la_SOURCES = \
                                SVN_Info.C \
                                Cut_Data.C \
                                Selector.C \
                                Variable_Selector.C \
                                Combined_Selector.C \
                                FastJet_Criterion.C \
                                Jet_Finder.C \
                                Jet_Selector.C \
                                Isolation_Selector.C \
                                NJet_Finder.C \
                                Standard_Selector.C \
                                Decay_Selector.C \
                                Decay2_Selector.C \
                                MinSelector.C \
                                DressedParticleSelector.C \
                                Fastjet_Finder.C \
                                Fastjet_Selector.C \
				Fastjet_Veto.C

localinc_HEADERS      = \
                        Cut_Data.H \
                        Selector.H \
                        Combined_Selector.H \
                        Jet_Finder.H

libPhasicSelectors_la_LIBADD = \
                @CONDITIONAL_FASTJETLIBS@

libPhasicSelectors_la_CPPFLAGS = $(AM_CPPFLAGS) \
                @CONDITIONAL_FASTJETINCS@
