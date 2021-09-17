#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#define COMPILE__Getter_Function
#define PARAMETER_TYPE CFPSHOWER::Kernel_Info
#define OBJECT_TYPE    CFPSHOWER::Kinematics_Base
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

Kinematics_Base::Kinematics_Base(const Kernel_Info & info) :
  m_name("unknown"), m_scheme(kin_type::none)
{}

Kinematics_Base::~Kinematics_Base() {}

void Kinematics_Base::InitSystem(Splitting & split, const Mass_Selector * msel)  {
  m_psplit  = split.GetSplitter()->Mom();
  m_pspect  = split.GetSpectator()->Mom();
  m_pboth   = m_psplit + m_pspect;
  m_msplit  = msel->Mass(split.GetSplitter()->Flav());
  m_msplit2 = sqr(m_msplit);
  m_mspect  = msel->Mass(split.GetSpectator()->Flav());
  m_mspect2 = sqr(m_mspect);
  m_Q2      = m_pboth.Abs2();
  m_Q       = sqrt(m_Q2);
  m_ismassive = (m_mspect>0. || m_m[0]>0. || m_m[1]>0.);
  m_allmomenta.clear();
  m_allmasses2.clear();
}

double Kinematics_Base::
Lambda(const double & a,const double & b,const double & c) {
  double lambda2 = Lambda2(a,b,c); 
  if (lambda2<0.) {
    msg_Error()<<"Error in "<<METHOD
	       <<"("<<a<<", "<<b<<", "<<c<<") yields nan.\n"
	       <<"   return 0. and hope for the best.\n";
    return 0.;
  }
  return sqrt(lambda2);
}

double Kinematics_Base::
Lambda2(const double & a,const double & b,const double & c) {
  return sqr(a-b-c)-4.*b*c;
}

