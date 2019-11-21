#include "MCATNLO/Showers/SF_Coupling.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE MCATNLO::SF_Key
#define OBJECT_TYPE MCATNLO::SF_Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

template class ATOOLS::Getter_Function
<MCATNLO::SF_Coupling,MCATNLO::SF_Key,SORT_CRITERION>;

template class ATOOLS::Getter_Function
<void,MCATNLO::SFC_Filler_Key,SORT_CRITERION>;

using namespace MCATNLO;

double SF_Coupling::s_qfac=1.0;

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

void SF_Coupling::ColorPoint(Parton *const p) const
{
}

double SF_Coupling::ColorWeight(const Color_Info &ci) const
{
  return 1.0;
}
