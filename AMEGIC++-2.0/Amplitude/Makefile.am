SUBDIRS = Zfunctions AmplTools
lib_LTLIBRARIES = libAmplitude.la
libAmplitude_la_SOURCES = Single_Amplitude.C Amplitude_Manipulator.C \
                          Amplitude_Handler.C Prop_Generator.C \
			  Amplitude_Generator.C Color_Generator.C \
	                  Zfunc_Generator.C Amplitude_Group.C \
			  Single_Amplitude_Base.C Super_Amplitude.C Zfunc.C	
INCLUDES = -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Main \
	   -I../Model \
	   -I../String \
	   -I../ISR \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = Single_Amplitude.H \
                 Amplitude_Generator.H Amplitude_Manipulator.H \
		 Color_Generator.H  Amplitude_Handler.H \
	         Zfunc_Generator.H Pfunc.H Prop_Generator.H \
		 Zfunc.H Amplitude_Base.H Amplitude_Group.H \
		 Color_Group.H Super_Amplitude.H Single_Amplitude_Base.H \
		 Super_Amplitude.H
