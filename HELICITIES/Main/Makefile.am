pkglib_LTLIBRARIES       =	libHelicitiesMain.la

libHelicitiesMain_la_SOURCES       =	\
			C_Tensor.C \
			C_Spinor.C \
			C_Vector.C \
			Polarization_Tools.C  \
			XYZFuncs.C  \
			Spin_Structure.C \
			Partial_Amplitude_Base.C \
			Full_Amplitude_Base.C \
			Three_Particle_Amplitudes.C \
			Four_Particle_Amplitudes.C

pkginclude_HEADERS    =	\
			C_Tensor.H \
			C_Spinor.H \
			C_Vector.H \
			Polarization_Tools.H  \
			XYZFuncs.H \
			Spin_Structure.H \
			Partial_Amplitude_Base.H \
			Full_Amplitude_Base.H \
			Three_Particle_Amplitudes.H \
			Four_Particle_Amplitudes.H


AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      = @HELICITIESINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
