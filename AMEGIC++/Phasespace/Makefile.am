pkglib_LTLIBRARIES = libAmegicPSGen.la

SVNTAG = AMEGIC++/Phasespace
include ../../svn.make

libAmegicPSGen_la_SOURCES = SVN_Info.C \
			    Phase_Space_Generator.C Channel_Generator_Base.C \
			    Channel_Generator.C Channel_Generator_NPV.C \
                            Channel_Generator3V.C Channel_Generator3_NPV.C \
	 		    Channel_Generator_Decays.C Channel_Generator_KK.C \
			    Channel_Generator_UniV.C

localinc_HEADERS   = Phase_Space_Generator.H Channel_Generator_Base.H \
		 Channel_Generator.H Channel_Generator_NPV.H \
		 Channel_Generator3V.H Channel_Generator3_NPV.H \
		 Channel_Generator_Decays.H Channel_Generator_KK.H \
		 Channel_Generator_UniV.H
