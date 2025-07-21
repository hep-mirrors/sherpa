#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

HR_Parameters::HR_Parameters()  {}

HR_Parameters::~HR_Parameters() {}

double HR_Parameters::operator()(const std::string& keyword) const {
  return 0.;
}

HR_Resonance::HR_Resonance(Flavour flav,Flavour fl1,Flavour fl2,
			   const double & symmfac, const double & BR,
			   const double & mass, const double & width) :
  m_flav(flav), m_BR(BR), m_m2(sqr(fl1.HadMass()+fl2.HadMass())), m_symmfac(symmfac), 
  m_mass(mass>0.?mass:m_flav.HadMass()), m_mass2(m_mass*m_mass), 
  m_width(width>0.?width:m_flav.Width()), m_width2(m_width*m_width),
  m_prest(sqrt(m_mass2-m_m2)),
  m_spin(int(m_flav.Spin())), m_spin1(int(fl1.Spin())), m_spin2(int(fl2.Spin())),
  m_spinfac(double(2*m_spin+1)/double((2*m_spin1+1)*(2*m_spin2+1)))
{}

const double HR_Resonance::XStot_R(const double & s) const {
  if (s<m_m2) return 0.;
  double width2 = Width2(s);
  double xstot = (4.*M_PI/(s-m_m2) * m_symmfac * m_spinfac *
		  s * m_BR * width2 / (sqr(s-m_mass2)+s*width2));
  msg_Out()<<METHOD<<"(R = "<<m_flav<<", E = "<<sqrt(s)<<"): "<<m_BR<<", "
	   <<width2<<" ("<<m_width2<<") from "
	   <<m_prest<<" --> "<<xstot<<".\n";
  return xstot;
}

const double HR_Resonance::Width2(const double & s) const {
  if (s<m_m2) return 0.; 
  return m_width2*pow(sqrt(s-m_m2)/m_prest,2*m_spin+1);
}

