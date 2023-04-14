#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/Main/Coulomb.H"


using namespace MODEL;
using namespace YFS;

Coulomb::Coulomb(){
  m_MW = Flavour(kf_Wplus).Mass();
  m_GW = Flavour(kf_Wplus).Width();
  m_cms = rpa->gen.Ecms();
  m_subtract = 0;
}


Coulomb::~Coulomb(){
  
}


void Coulomb::Calculate(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2){
  m_s1 = p1.Abs2();
  m_s2 = p2.Abs2();
  // PRINT_VAR(m_s);
  //Eq 9 in https://arxiv.org/pdf/hep-ph/9507422.pdf
  double E = (m_s-4*sqr(m_MW))/(4*m_MW);
  double sarg = sqrt(sqr(E)+sqr(m_GW));
  double pp1 = sqrt(0.5*m_MW*(sarg-E));
  double pp2 = sqrt(0.5*m_MW*(sarg+E));
  double pp = 1./(4.*m_s)*(sqr(m_s) -2*m_s*(m_s1+m_s2) +sqr(m_s1-m_s2));
  double p  = sqrt(pp);
  m_p = p;
  double absKappa = m_MW*sqrt(sqr(E)+sqr(m_GW));
  double arg = (absKappa-pp)/(2*p*pp1);
  m_weight = 1 + m_alpha*sqrt(m_s)/(4*p)*(M_PI - 2.*atan(arg));
  if(IsBad(m_weight)){
    msg_Error()<<METHOD<<std::endl
                   <<"E = "<< E<<std::endl
                   <<"sqrt(s) = "<< sqrt(m_s)<<std::endl
                   <<"sqrt(s1) = "<< sqrt(m_s1)<<std::endl
                   <<"sqrt(s2) = "<< sqrt(m_s1-m_s2)<<std::endl
                   <<"p1 = "<< pp1<<std::endl
                   <<"p2 = "<< pp2<<std::endl
                   <<"pp = "<< pp<<std::endl
                   <<"k = "<< absKappa<<std::endl
                   <<"arg = "<< arg<<std::endl
                   <<"arctan(arg) = "<< atan(arg)<<std::endl
                   <<"weight = "<< m_weight<<std::endl;
  }
}


void Coulomb::Subtract(){
  m_weight-= sqrt(m_s)*m_alpha/(2*m_p)*M_PI;
  // m_weight-= 0;
}
