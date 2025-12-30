#include "BEAM/Spectra/Fixed_Target.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

#include <string>
#include <fstream>

using namespace BEAM;
using namespace ATOOLS;


Fixed_Target::Fixed_Target(const ATOOLS::Flavour _beam,const double &_energy,
	    const double &_polarisation, const int &_fixed) :
  Beam_Base(beamspectrum::Fixed_Target,_beam,_energy,_polarisation,_fixed)
  {
  	if(_fixed == 1) m_beam_lab = m_lab;
  	else {
  		m_fixed_target = Vec4D(_beam.Mass(),0,0,0);
  		m_lab = m_fixed_target;
  	}
  	m_on = false;
  }

Fixed_Target::~Fixed_Target(){}

bool Fixed_Target::CalculateWeight(double x,double q2)  {return 1;}
// double Fixed_Target::Weight(Flavour fl)                { return m_weight; }
// ATOOLS::Flavour Fixed_Target::Remnant()                { return kf_none; }


