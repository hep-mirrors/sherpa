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

std::ostream & CFPSHOWER::operator<<(std::ostream &s,const ffff_mode::code & mode) {
  s<<(mode==ffff_mode::same?"same":"diff");
  return s;
}

std::ostream & CFPSHOWER::operator<<(std::ostream &s,const subtract::code & sub) {
  if      (sub==subtract::none) s<<"none";
  else if (sub==subtract::coll) s<<"coll";
  else if (sub==subtract::soft) s<<"soft";
  else if (sub==subtract::both) s<<"both";
  return s;
}

double SF_Base::Nc = 3.;
double SF_Base::CF = (ATOOLS::sqr(SF_Base::Nc)-1.)/(2.*SF_Base::Nc);
double SF_Base::CA = SF_Base::Nc;
double SF_Base::TR = 1./2.;

SF_Base::SF_Base(const Kernel_Info & info) :
  m_split(info.GetSplit()),
  m_flavs(info.GetFlavs()), m_tagsequence(info.TagSequence()),
  m_name("generic SF"),
  m_CMW(info.KFactor())
{ }  

double SF_Base::Lambda(const double & a,const double & b,const double & c) const {
  double lambda2 = Lambda2(a,b,c); 
  if (lambda2<0.) {
    msg_Error()<<"Error in "<<METHOD<<"("<<a<<", "<<b<<", "<<c<<") yields nan.\n"
	       <<"   return 0. and hope for the best.\n";
    return 0.;
  }
  return sqrt(lambda2);
}

double SF_Base::Lambda2(const double & a,const double & b,const double & c) const {
  return sqr(a-b-c)-4.*b*c;
}
