#include "BEAM/Main/Monochromatic.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace BEAM;
using namespace std;

Monochromatic::Monochromatic(const Flavour _beam,const double _energy,
			     const double _polarisation,const int _dir) :
  Beam_Base(string("Monochromatic"),_beam,_energy,_polarisation,_dir)
{ }


Beam_Base * Monochromatic::Copy() 
{
  return new Monochromatic(m_beam,m_energy,m_polarisation,m_dir);
}

bool Monochromatic::CalculateWeight(double x,double q2) { return 1; }
double Monochromatic::Weight(Flavour fl)                { return m_weight; }
ATOOLS::Flavour Monochromatic::Remnant()                { return kf_photon; }





