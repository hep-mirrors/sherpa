pkglib_LTLIBRARIES       =	libPDF.la

libPDF_la_SOURCES     =	\
		    	PDF_Base.C \
		    	PDF_Electron.C \
			Continued_PDF.C \
		    	PDF_Handler.C \
			ISR_Base.C \
		    	Intact.C \
			Structure_Function.C \
			ISR_Handler.C 

pkginclude_HEADERS       =	\
	         	PDF_Base.H \
		 	PDF_Electron.H \
			Continued_PDF.H \
			PDF_Handler.H \
			ISR_Base.H \
			Intact.H \
			Structure_Function.H\
			ISR_Info.H \
		 	ISR_Handler.H 

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@SHERPAINCS@ @APACICINCS@ @PDFINCS@ @CTEQINCS@ \
		@MODELINCS@ @BEAMINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
