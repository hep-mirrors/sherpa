#include "BEAM/Spectra/DM_beam.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

#include <string>
#include <fstream>

using namespace BEAM;
using namespace ATOOLS;

std::ostream & BEAM::operator<<(std::ostream &s,const DM_type::code & type) {
  switch (int(type)) {
  case 1: s<<"Boltzmann";break;
  case 0:
  default: s<<"Unknown";break;
  }
  return s;
}

DM_beam::DM_beam(const Flavour beam,const double & temperature,
		 const DM_type::code & form, const bool & relativistic,
		 const int dir):
  Beam_Base(beamspectrum::DM,beam,1.e6,0.,dir),
  m_mass(beam.Mass(true)), m_temperature(temperature),
  m_formfactor(form), m_relativistic(relativistic), m_debugging(false)
{
  if (m_debugging) {
    std::string filename = "DM_"+ToString(int((1+dir)/2))+".log";
    selfTest(filename);
  }
}

bool DM_beam::CalculateWeight(double x,double s)
{
  // event generation mode
  double E = x*sqrt(s);
  if (E < m_mass) {
    m_weight = 0;
    return true;
  }
  double p = sqrt(sqr(E) - sqr(m_mass));
  m_weight = E*p * exp(-(E-m_mass)/m_temperature);
  return true;
}

void DM_beam::selfTest(std::string filename) {}
