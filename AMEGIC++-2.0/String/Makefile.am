lib_LTLIBRARIES = libString.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

libString_la_SOURCES = String_Handler.C String_Generator.C \
                       String_Tree.C String_Output.C \
	               MyString.C String_Library.C
INCLUDES = -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Main \
	   -I../Model \
	   -I../String \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
	   -I../../MODEL-1.0/Main 
noinst_HEADERS = String_Handler.H String_Generator.H \
                 String_Tree.H String_Output.H \
	         MyString.H Kabbala.H Values.H String_Library.H

EXTRA_DIST    = $(CVSINCLUDE)