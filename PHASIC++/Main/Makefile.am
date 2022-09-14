
pkglib_LTLIBRARIES    =	libPhasicMain.la

GITTAG = PHASIC++/Main
include ../../git.make

libPhasicMain_la_SOURCES      =	\
				Git_Info.C \
				Color_Integrator.C \
				Helicity_Integrator.C \
				Process_Integrator.C \
				Phase_Space_Handler.C \
				Phase_Space_Integrator.C \
				Event_Reader.C

localinc_HEADERS      =	\
			Color_Integrator.H \
			Helicity_Integrator.H \
		 	Process_Integrator.H \
			Phase_Space_Handler.H \
			Phase_Space_Integrator.H \
			Event_Reader.H

dist_bin_SCRIPTS = plot_stats.sh
