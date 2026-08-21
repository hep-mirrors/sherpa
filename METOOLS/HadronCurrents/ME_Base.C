#include "METOOLS/HadronCurrents/ME_Base.H"
#include "ATOOLS/Org/Message.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE    METOOLS::ME_Base
#define PARAMETER_TYPE METOOLS::ME_Parameters
#include "ATOOLS/Org/Getter_Function.C"

using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

ME_Base::ME_Base(const ATOOLS::Flavour_Vector& flavs,
		 const std::vector<int>& decayindices,
		 const std::string& name) :
  METOOLS::Spin_Amplitudes(flavs,Complex(0.0,0.0)),
  m_name(name), m_flavs(flavs)
{

  m_pi.resize(m_flavs.size());
  m_masses.resize(m_flavs.size());
  m_masses2.resize(m_flavs.size());
  for(int meindex=0; meindex<m_flavs.size(); meindex++) {
    m_pi[meindex]      = decayindices[meindex];
    m_masses[meindex]  = m_flavs[m_pi[meindex]].HadMass();
    m_masses2[meindex] = sqr(m_masses[meindex]);
  }
}

double ME_Base::lambdaNorm(const double M,const double m1,const double m2) {
  return sqrt((sqr(M)-sqr(m1+m2))*(sqr(M)-sqr(m1-m2)))/(2.*M);
}

void ME_Base::SetModelParameters(ATOOLS::Scoped_Settings& s)
{
  DEBUG_FUNC("");
  GeneralModel model;
  for (const auto& key: s.GetKeys()) {
    //if (m_name=="Current_ME" && key=="J1") continue;
    //if (m_name=="Current_ME" && key=="J2") continue;
    if (s[key].GetItemsCount()>1) {
      model[key+string("_abs")] = s[key][0].GetScalarWithOtherDefault(-1.0);
      model[key+string("_phase")] = s[key][1].GetScalarWithOtherDefault(0.0);
    }
    else model[key] = s[key].GetScalarWithOtherDefault(-1.0);
  }
  SetModelParameters(model);
}
