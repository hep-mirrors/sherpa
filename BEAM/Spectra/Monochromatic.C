#include "BEAM/Spectra/Monochromatic.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace BEAM;
using namespace std;

Monochromatic::Monochromatic(const Flavour _beam,const double _energy,
			     const double _polarisation,const int _dir) :
  Beam_Base(beamspectrum::monochromatic,_beam,_energy,_polarisation,_dir)
{ }


Beam_Base * Monochromatic::Copy()
{
  return new Monochromatic(m_beam,m_energy,m_polarisation,m_dir);
}

ATOOLS::Flavour Monochromatic::Remnant()                { return kf_none; }





