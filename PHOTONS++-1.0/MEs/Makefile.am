pkglib_LTLIBRARIES = libPhotonsMEs.la

libPhotonsMEs_la_SOURCES = 	PHOTONS_ME_Base.C \
				Vector_To_Fermion_Fermion.C \
				Scalar_To_Fermion_Fermion.C \
				Scalar_To_Scalar_Scalar.C \
				Z_To_Lepton_Lepton.C \
				W_To_Lepton_Neutrino.C \
				Tau_To_Lepton_Neutrinos.C

pkginclude_HEADERS = 	ME_List.H \
			PHOTONS_ME_Base.H \
			Vector_To_Fermion_Fermion.H \
			Scalar_To_Fermion_Fermion.H \
			Scalar_To_Scalar_Scalar.H \
			Z_To_Lepton_Lepton.H \
			W_To_Lepton_Neutrino.H \
			Tau_To_Lepton_Neutrinos.H

AM_CXXFLAGS = @SHERPAFLAGS@

INCLUDES = @ATOOLSINCS@ @PHOTONSINCS@ @MODELINCS@ @HELICITIESINCS@

EXTRA_DIST    = @SVNINCLUDE@
