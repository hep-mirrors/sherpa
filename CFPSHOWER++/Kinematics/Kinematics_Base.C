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

double Kinematics_Base::Lambda(const double & a,const double & b,const double & c) {
  double lambda2 = Lambda2(a,b,c); 
  if (lambda2<0.) {
    msg_Error()<<"Error in "<<METHOD<<"("<<a<<", "<<b<<", "<<c<<") yields nan.\n"
	       <<"   return 0. and hope for the best.\n";
    return 0.;
  }
  return sqrt(lambda2);
}

double Kinematics_Base::Lambda2(const double & a,const double & b,const double & c) {
  return sqr(a-b-c)-4.*b*c;
}

