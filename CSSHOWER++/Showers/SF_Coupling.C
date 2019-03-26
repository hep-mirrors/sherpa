#include "CSSHOWER++/Showers/SF_Coupling.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE CSSHOWER::SF_Key
#define OBJECT_TYPE CSSHOWER::SF_Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

template class ATOOLS::Getter_Function
<CSSHOWER::SF_Coupling,CSSHOWER::SF_Key,SORT_CRITERION>;

template class ATOOLS::Getter_Function
<void,CSSHOWER::SFC_Filler_Key,SORT_CRITERION>;

using namespace CSSHOWER;

SF_Coupling::SF_Coupling(const SF_Key &key):
  p_lf(NULL), m_type(key.m_type),
  m_cplfac(1.0), m_kfmode(key.m_kfmode)
{
}

SF_Coupling::~SF_Coupling() {}

double SF_Coupling::CplFac(const double &scale) const
{
  return m_cplfac;
}
