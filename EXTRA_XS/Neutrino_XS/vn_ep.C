#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"


using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;
using namespace PHASIC;
using namespace std;

namespace EXTRAXS {

  class vn_ep : public ME2_Base
  {
  private:
  public:
    vn_ep(const External_ME_Args &args, const incomingboson::code &boson, const incomingnucleon::code &nucleon);

    double operator()(const ATOOLS::Vec4D_Vector &mom);
    double m_alpha, m_alphas, m_s;
    Flavour m_flv;
    std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
    incomingboson::code m_boson;
    incomingnucleon::code m_nucleon;
  };

  vn_ep::vn_ep(const External_ME_Args &args, const incomingboson::code &boson, const incomingnucleon::code &nucleon)
      : ME2_Base(args), m_boson(boson), m_nucleon(nucleon)
  {
    m_oew = 0;
    m_oqcd = 0;
    m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");
    m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");

    p_formfactor = std::unique_ptr<FormFactor_EMnucleon>(new FormFactor_EMnucleon(m_boson, m_nucleon));
  }

  double vn_ep::operator()(const ATOOLS::Vec4D_Vector &momenta)
  {
    
    double m1 = momenta[0].Mass();
    m_s = (momenta[2] + momenta[3]).Abs2();
    double t((momenta[0] - momenta[2]).Abs2());
    double amp = 1;
    return 32 * sqr((*aqed)(m_s)) * amp / m_s / m_s; //*p_formfactor->Eval(s); // flux = 0.5
  }

}  // namespace EXTRAXS

DECLARE_TREEME2_GETTER(EXTRAXS::vn_ep, "vn_ep")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::vn_ep>::operator()(const External_ME_Args &args) const
{
  return new vn_ep(args, incomingboson::W, incomingnucleon::neutron);
}
