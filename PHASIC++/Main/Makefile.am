
pkglib_LTLIBRARIES    =	libPhasicMain.la

GITTAG = PHASIC++/Main
include ../../git.make

libPhasicMain_la_SOURCES      =	\
				Git_Info.C \
				Color_Integrator.C \
				Helicity_Integrator.C \
				Process_Integrator.C \
				Channel_Creator.C \
				Phase_Space_Point.C \
				Phase_Space_Enhance.C \
				Phase_Space_Handler.C \
				Phase_Space_Integrator.C

localinc_HEADERS      =	\
			Color_Integrator.H \
			Helicity_Integrator.H \
		 	Process_Integrator.H \
			Channel_Creator.H \
			Phase_Space_Point.H \
			Phase_Space_Enhance.H \
			Phase_Space_Variation_Weights.H \
			Phase_Space_Handler.H \
			Phase_Space_Integrator.H

dist_bin_SCRIPTS = plot_stats.sh
