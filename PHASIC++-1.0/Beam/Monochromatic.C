#include "Monochromatic.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace BEAM;
using namespace std;

Monochromatic::Monochromatic(const Flavour _beam,const double _energy,
			     const double _polarization,bool & okay) 
{
  m_beam         = _beam;
  m_energy       = dabs(_energy);
  m_polarization = _polarization;
  m_type         = string("Monochromatic");
  m_weight       = 1.;

  double disc    =  1.-sqr(m_beam.PSMass()/m_energy);
  if (disc<0) {
    msg.Error()<<"Error in Monochromatic::Monochromatic :"<<endl
	       <<"   Mismatch of energy and mass of beam particle : "<<m_beam<<" / "<<m_energy<<endl
	       <<"   Will lead to termination of program."<<endl;
    okay         = 0;
    return;
  }
  m_lab          = Vec4D(m_energy,0.,0.,_energy*sqrt(disc));
  okay           = 1;
  msg.Tracking()<<"Initialised Monochromatic for beam "<<m_beam<<endl;
}

Beam_Base * Monochromatic::Copy() 
{
  double energy = Sign(int(m_lab[3]))*m_energy;
  bool okay     = 1;
  return new Monochromatic(m_beam,energy,m_polarization,okay);
}

bool Monochromatic::CalculateWeight(double x,double q2) { return 1; }
double Monochromatic::Weight(Flavour fl)                { return m_weight; }





