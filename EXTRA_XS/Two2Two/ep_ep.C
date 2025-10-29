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

    class ep_ep : public ME2_Base
    {
    private:
    public:
        ep_ep(const External_ME_Args &args, const incomingboson::code &boson, const incomingnucleon::code &nucleon);

        double operator()(const ATOOLS::Vec4D_Vector &mom);
        double m_alpha, m_alphas, m_s;
        Flavour m_flv;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
        incomingboson::code m_boson;
        incomingnucleon::code m_nucleon;
    };

    ep_ep::ep_ep(const External_ME_Args &args, const incomingboson::code &boson, const incomingnucleon::code &nucleon)
        : ME2_Base(args), m_boson(boson), m_nucleon(nucleon)
    {
        m_oew = 0;
        m_oqcd = 0;
        m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");
        m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");

        p_formfactor = std::unique_ptr<FormFactor_EMnucleon>(new FormFactor_EMnucleon(m_boson, m_nucleon));
    }

    double ep_ep::operator()(const ATOOLS::Vec4D_Vector &momenta)
    {
        //double m1 = momenta[0].Mass();
        //m_s = (momenta[2] + momenta[3]).Abs2();
        //double t((momenta[0] - momenta[2]).Abs2());
        double Q2 = (momenta[0] - momenta[2]).Abs2();
        NucleonFormFactors ff = p_formfactor->GetFormFactors(Q2);
        double F1 = ff.F1;
        double F2 = ff.F2;
        double amp = F1 * F2;
        return 32 * amp; //*p_formfactor->Eval(s); // flux = 0.5
    }
}

DECLARE_TREEME2_GETTER(EXTRAXS::ep_ep, "ep_ep")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::ep_ep>::
operator()(const External_ME_Args &args) const
{
    const Flavour_Vector fl = args.Flavours();
    if (fl.size() != 4) return NULL;
    
    // Check for correct process e- p -> e- p 
    if (fl[0] == Flavour(kf_e) && fl[1] == Flavour(kf_p_plus) &&
        fl[2] == Flavour(kf_e) && fl[3] == Flavour(kf_p_plus))
    {
        return new ep_ep(args, incomingboson::photon, incomingnucleon::proton);
    }
    
    return NULL;
}
