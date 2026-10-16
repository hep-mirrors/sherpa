#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "METOOLS/HadronCurrents/FormFactors/Parameter_Sets.H"

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
  CheckReservedModel(m_ffmodel);
  // Two loops, NOT one: m_pi holds the decay indices and is generally SHORTER
  // than m_flavs (five flavours but three hadrons for e+e- -> pi+ pi- pi0, or
  // for a tau decay).  Indexing m_pi by the flavour loop reads past its end --
  // harmless for a current that only ever touches m_pi[0..2], but it fills the
  // vector with garbage that bites anything iterating over m_pi.size().
  for (size_t i=0;i<params.m_flavs.size();i++) {
    m_flavs.push_back(params.m_flavs[i]);
    m_masses.push_back(params.m_masses[i]);
    m_masses2.push_back(params.m_masses2[i]);
  }
  for (size_t i=0;i<params.m_pi.size();i++) {
    m_pi.push_back(params.m_pi[i]);
  }

  // Local parameter sets; see Parameter_Sets.H.  Applied HERE rather than in
  // the caller so that every consumer of a form factor gets them, decay or
  // production alike.  ff_model::none means "no dynamical model" -- a constant
  // form factor, or a current whose card named no FORM_FACTOR -- which no
  // published fit can describe, so those are left alone.
  if (m_ffmodel!=ff_zero && m_ffmodel!=ff_one)
    ApplyParameterSet(p_model,m_flavs,m_pi,m_ffmodel);
}


void FormFactor_Base::CheckReservedModel(int model) {
  if (model==ff_invalid || model<0 || (model>1 && model<10))
    THROW(fatal_error,"Invalid FORM_FACTOR model number " +
          std::to_string(model) + ". Use 0 (zero), 1 (one), or a family model >= 10; 99 is reserved.");
}
