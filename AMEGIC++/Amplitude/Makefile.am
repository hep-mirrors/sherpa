SUBDIRS = Zfunctions

pkglib_LTLIBRARIES = libAmplitude.la

SVNTAG = AMEGIC++/Amplitude
include ../../svn.make

libAmplitude_la_SOURCES = SVN_Info.C \
			  Single_Vertex.C Vertex.C \
			  CFColor.C Amplitude_Base.C \
			  Single_Amplitude.C Amplitude_Manipulator.C \
                          Amplitude_Generator.C Color_Generator.C \
	                  Prop_Generator.C Zfunc_Generator.C Amplitude_Group.C \
		          Amplitude_Handler.C \
			  Single_Amplitude_Base.C Super_Amplitude.C Zfunc.C \
                          FullAmplitude_MHV_Base.C \
                          FullAmplitude_External.C \
                          Amplitude_Output.C Lorentz_Function.C Lorentz_Functions.C

localinc_HEADERS   = \
		 Single_Vertex.H Vertex.H \
		 CFColor.H Single_Amplitude.H \
                 Amplitude_Generator.H Amplitude_Manipulator.H \
		 Color_Generator.H \
	         Prop_Generator.H Zfunc_Generator.H Pfunc.H \
		 Zfunc.H Amplitude_Handler.H Amplitude_Base.H \
		 Amplitude_Group.H \
		 Color_Group.H Super_Amplitude.H Single_Amplitude_Base.H \
		 FullAmplitude_MHV_Base.H \
		 FullAmplitude_External.H \
		 Amplitude_Output.H Lorentz_Function.H
