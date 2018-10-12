bin_PROGRAMS  =	Sherpa@VERSIONING@

Sherpa@VERSIONING@_SOURCES        =	Main.C 

Sherpa@VERSIONING@_LDADD  =	@SHERPABUILDDIR@/Main/libSherpaMain.la @ATOOLSLIBS@ \
                @CONDITIONAL_FASTJETLIBS@

dist_pkgdata_DATA         = \
			Analysis.dat \
			Fragmentation.dat \
			Integration.dat \
			ISR.dat \
			ME.dat \
			MI.dat \
			Run.dat \
			Shower.dat

bin_SCRIPTS = Sherpa@VERSIONING@-config
dist_pkgdata_SCRIPTS 	= \
			get_bibtex \
			plot_graphs
