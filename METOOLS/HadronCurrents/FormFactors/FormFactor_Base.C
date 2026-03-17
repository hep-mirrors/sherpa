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
  m_flavs(params.m_flavs),
  m_masses(params.m_masses), m_masses2(params.m_masses2), m_pi(params.m_pi),
  m_name(params.m_name),
  p_model(params.p_model)
{ }

