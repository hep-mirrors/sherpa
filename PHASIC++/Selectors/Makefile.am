
pkglib_LTLIBRARIES    =	libPhasicSelectors.la

SVNTAG = PHASIC++__Selectors
include ../../svn.make

libPhasicSelectors_la_SOURCES =	\
				SVN_Info.C \
				Cut_Data.C \
				Selector.C \
				Variable_Selector.C \
				Combined_Selector.C \
				Jet_Finder.C \
				KT_Finder.C \
				NJet_Finder.C \
				Cone_Finder.C \
				Standard_Selector.C \
				NLO_Selector.C \
				MinSelector.C \
				Selector_Bias.C \
				Fastjet_Finder.C

localinc_HEADERS      =	\
			SVN_Info.C.in \
			Cut_Data.H \
			Selector.H \
			Combined_Selector.H \
			Jet_Finder.H \
			KT_Finder.H

libPhasicSelectors_la_LIBADD = \
		@CONDITIONAL_FASTJETLIBS@
        
INCLUDES      = @CONDITIONAL_FASTJETINCS@
