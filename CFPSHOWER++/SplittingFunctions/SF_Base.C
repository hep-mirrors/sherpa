#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#define COMPILE__Getter_Function
#define PARAMETER_TYPE CFPSHOWER::Kernel_Info
#define OBJECT_TYPE    CFPSHOWER::SF_Base
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

std::ostream & CFPSHOWER::operator<<(std::ostream &s,const subtract::code & sub) {
  if      (sub==subtract::none) s<<"none";
  else if (sub==subtract::coll) s<<"coll";
  else if (sub==subtract::soft) s<<"soft";
  else if (sub==subtract::both) s<<"both";
  return s;
}

SF_Base::SF_Base(const Kernel_Info & info) :
  m_split(info.GetSplit()), 
  m_flavs(info.GetFlavs()), m_tags(info.TagSequence()), m_nout(m_flavs.size()),
  m_type(info.Type()),
  m_name("generic SF"),
  m_CMW(info.KFactor()), m_softcorr(info.SoftCorrection()),
  m_endpoint(info.Endpoint()), m_subtract(subtract::none),
  m_weight(0.)
{ }  
