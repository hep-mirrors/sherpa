#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE METOOLS::FormFactor_Base
#define PARAMETER_TYPE METOOLS::FF_Parameters
#include "ATOOLS/Org/Getter_Function.C"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

FormFactor_Base::FormFactor_Base(const FF_Parameters & params) :
  m_ffmodel(params.m_ffmodel),
  m_name(params.m_name),
  p_model(params.p_model)
{
  for (size_t i=0;i<params.m_flavs.size();i++) {
    m_flavs.push_back(params.m_flavs[i]);
    m_masses.push_back(params.m_masses[i]);
    m_masses2.push_back(params.m_masses2[i]);
    m_pi.push_back(params.m_pi[i]);
  }
}

