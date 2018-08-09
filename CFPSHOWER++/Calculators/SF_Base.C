#include "CFPSHOWER++/Calculators/SF_Base.H"
#define COMPILE__Getter_Function
#define PARAMETER_TYPE CFPSHOWER::Kernel_Info
#define OBJECT_TYPE    CFPSHOWER::SF_Base
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

SF_Base::SF_Base(const Kernel_Info & info) :
  m_flavs(info.GetFlavs()),
  m_CF(4./3.), m_CA(3.), m_TR(1./2.), m_zeta3(1.202056903159594),
  m_orderA(1), m_orderB(1),
  m_swap(info.Swapped()), 
  m_name("generic SF")
{
  // KFactor: 1 = include A2
  // KFactor: 2 = include B2
  // KFactor: 4 = include A3
  switch (info.KFactor()) {
  case 7: m_orderA = 3; m_orderB = 2; break;
  case 5: m_orderA = 3; m_orderB = 1; break;
  case 3: m_orderA = 2; m_orderB = 2; break;
  case 1: m_orderA = 2; m_orderB = 1; break;
  case 0: m_orderA = 1; m_orderB = 1; break;
  default:
    msg_Error()<<"Error in "<<METHOD<<": invalid higher order setting, will use LO.\n";
    break;
  }
}  

double SF_Base::S1(const double & q2) const { return 0.; }

double SF_Base::S2(const double & q2) const { return 0.; }

double SF_Base::Lambda(const double & a,const double & b,const double & c) const {
  if (sqr(a-b-c)<4.*b*c) {
    msg_Error()<<"Error in "<<METHOD<<"("<<a<<", "<<b<<", "<<c<<") yields nan.\n"
	       <<"   return 0. and hope for the best.\n";
    return 0.;
  }
  return sqrt(sqr(a-b-c)-4.*b*c);
}

double SF_Base::Lambda2(const double & a,const double & b,const double & c) const {
  return sqr(a-b-c)-4.*b*c;
}
