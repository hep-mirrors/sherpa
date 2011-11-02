pkglib_LTLIBRARIES = libPhotonsMEs.la

SVNTAG = PHOTONS++/MEs
include ../../svn.make

libPhotonsMEs_la_SOURCES =	SVN_Info.C \
				PHOTONS_ME_Base.C \
				Scalar_To_Scalar_Scalar.C \
				Scalar_To_Fermion_Fermion.C \
				Scalar_To_Scalar_Lepton_Neutrino.C \
				Scalar_To_Vector_Lepton_Neutrino.C \
				Tau_To_Lepton_Neutrinos.C \
				W_To_Lepton_Neutrino.C \
				Z_To_Fermion_Fermion.C \
				Vector_To_Fermion_Fermion.C \
				Vector_To_Scalar_Scalar.C

libPhotonsMEs_la_LIBADD = @METOOLSLIBS@

localinc_HEADERS =	PHOTONS_ME_Base.H \
			Scalar_To_Scalar_Scalar.H \
			Scalar_To_Fermion_Fermion.H \
			Scalar_To_Scalar_Lepton_Neutrino.H \
			Scalar_To_Vector_Lepton_Neutrino.H \
			Tau_To_Lepton_Neutrinos.H \
			W_To_Lepton_Neutrino.H \
			Z_To_Fermion_Fermion.H \
			Vector_To_Fermion_Fermion.H \
			Vector_To_Scalar_Scalar.H
