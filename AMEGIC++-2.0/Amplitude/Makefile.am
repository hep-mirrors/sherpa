SUBDIRS = Zfunctions
pkglib_LTLIBRARIES = libAmplitude.la

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

AMEGICFLAGS   =	# -pedantic -Wall -Winline
DEFS          =

libAmplitude_la_SOURCES = CFColor.C Amplitude_Base.C \
			  Single_Amplitude.C Amplitude_Manipulator.C \
                          Amplitude_Generator.C Color_Generator.C \
	                  Prop_Generator.C Zfunc_Generator.C Amplitude_Group.C \
		          Amplitude_Handler.C \
			  Single_Amplitude_Base.C Super_Amplitude.C Zfunc.C \
                          Amplitude_Output.C
INCLUDES = -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Main \
	   -I../Model \
	   -I../String \
	   -I../ISR \
	   -I../../MODEL-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
pkginclude_HEADERS = CFColor.H Single_Amplitude.H \
                 Amplitude_Generator.H Amplitude_Manipulator.H \
		 Color_Generator.H \
	         Prop_Generator.H Zfunc_Generator.H Pfunc.H \
		 Zfunc.H Amplitude_Handler.H Amplitude_Base.H \
		 Amplitude_Group.H \
		 Color_Group.H Super_Amplitude.H Single_Amplitude_Base.H \
		 Super_Amplitude.H Amplitude_Output.H

EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(AMEGICFLAGS)
