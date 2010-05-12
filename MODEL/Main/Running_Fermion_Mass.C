#include "MODEL/Main/Running_Fermion_Mass.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace MODEL;
using namespace ATOOLS;

Running_Fermion_Mass::Running_Fermion_Mass(ATOOLS::Flavour _flav,double _polemass,
					   Running_AlphaS * _as) :
  m_polemass(_polemass), p_as(as)
{
  m_type    = std::string("Running Mass");
  m_name    = "Mass_"+ToString(_flav);
  m_defval  = m_polemass;
  if (_flav.Mass(true)<1.||(!_flav.IsQuark())||m_polemass<1.) {
    m_order = 0;
    p_as    = NULL;
    return;
  }
  double scale = sqr(m_polemass);
  m_order      = 1;
  m_gamma      = 1./M_PI;
  m_beta       = (33. - 2.*p_as->Nf(scale/4.))/(12.*M_PI); 
  m_alphahat   = (*p_as)(scale);
}

double Running_Fermion_Mass::operator()(double t) {
  if (m_order==0) return m_polemass;
  if (t<0.) t = -t;
  if (t<sqr(m_polemass)) return m_polemass;
  m_beta = (33. - 2.*p_as->Nf(t/4.))/(12.*M_PI); 

  return m_polemass*pow((*p_as)(t)/m_alphahat,m_gamma/m_beta);
}

void Running_Fermion_Mass::SelfTest() {
  double m_test = m_polemass/2.;
  for (int i=0;i<100;i++) {
    m_test += m_polemass/20.*i;
    std::cout<<"  "<<m_test<<" "<<(*this)(sqr(m_test))<<std::endl;
  }
}
