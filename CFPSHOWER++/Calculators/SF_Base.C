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

SF_Base::SF_Base(const Kernel_Info & info) :
  m_split(info.GetSplit()),
  m_flavs(info.GetFlavs()), m_tags(info.TagSequence()),
  m_name("generic SF"),
  m_CMW(info.KFactor()), m_subtract(subtract::none)
{
  m_invtags.resize(m_tags.size());
  for (size_t i=0;i<m_tags.size();i++) m_invtags[m_tags[i]] = i;
}  

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

