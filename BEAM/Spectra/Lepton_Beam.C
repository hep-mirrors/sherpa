#include "BEAM/Spectra/Lepton_Beam.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace BEAM;

Lepton_Beam::Lepton_Beam(const Flavour _beam,const double _energy,
			     const double _polarisation,const int _dir)
    : Beam_Base(beamspectrum::Leptonic,_beam,_energy,_polarisation,_dir),
      m_none(kf_none) {
      Settings &s = Settings::GetMainSettings();
      m_dev = s["BES"]["Deviation"].SetDefault(0.1).Get<double>();
      }


Beam_Base * Lepton_Beam::Copy()
{
  return new Lepton_Beam(m_beam,m_energy,m_polarisation,m_dir);
}

ATOOLS::Flavour Lepton_Beam::Remnant()                { return m_none; }





