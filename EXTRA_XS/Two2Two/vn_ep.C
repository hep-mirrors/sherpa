#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"
#include <memory>

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;
using namespace PHASIC;
using namespace std;

namespace EXTRAXS
{

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
    msg_Out() << "vn_ep::vn_ep(): Constructor called" << std::endl;
    m_oew = 0;
    m_oqcd = 0;
    msg_Out() << "vn_ep::vn_ep(): Getting alpha_QED" << std::endl;
    m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");
    msg_Out() << "vn_ep::vn_ep(): Getting strong_cpl" << std::endl;
    m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");

    msg_Out() << "vn_ep::vn_ep(): Creating FormFactor_EMnucleon" << std::endl;
    p_formfactor = std::make_unique<FormFactor_EMnucleon>(m_boson, m_nucleon);
    msg_Out() << "vn_ep::vn_ep(): Constructor complete" << std::endl;
  }

  double vn_ep::operator()(const ATOOLS::Vec4D_Vector &momenta)
  {
    double Q2 = (momenta[0] - momenta[2]).Abs2();
    NucleonFormFactors ff = p_formfactor->GetFormFactors(Q2);
    double F1 = ff.F1;
    double F2 = ff.F2;
    double amp = F1 * F2;
    return 32 * amp; //*p_formfactor->Eval(s); // flux = 0.5
  }
}

DECLARE_TREEME2_GETTER(EXTRAXS::vn_ep, "vn_ep")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::vn_ep>::
operator()(const External_ME_Args &args) const
{
  const Flavour_Vector fl = args.Flavours();
  if (fl.size() != 4)
    return NULL;

  // Check for ve n -> P+ e- 
  // Sherpas initial state: fl[0]=ve, fl[1]=n
  // Sherpas final state: fl[2]=P+, fl[3]=e-
  if (fl[0] == Flavour(kf_nue) && fl[1] == Flavour(kf_n) &&
      fl[2] == Flavour(kf_p_plus) && fl[3] == Flavour(kf_e))
  {
    return new vn_ep(args, incomingboson::W, incomingnucleon::neutron);
  }

  return NULL;
}