pkglib_LTLIBRARIES = libString.la

libString_la_SOURCES = String_Handler.C String_Generator.C \
                       String_Tree.C String_Output.C \
                       MyString.C Kabbala.C String_Library.C

pkginclude_HEADERS = String_Handler.H String_Generator.H \
                 String_Tree.H String_Output.H \
                 MyString.H Kabbala.H Values.H String_Library.H

INCLUDES = @AMEGICINCS@ @ATOOLSINCS@ @MODELINCS@

EXTRA_DIST    = @SVNINCLUDE@
AM_CXXFLAGS   = @SHERPAFLAGS@
