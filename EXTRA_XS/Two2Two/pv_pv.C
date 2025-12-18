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
    // NOTE: that this class also does nv-nv scattering when the nucleon is a neutron
    //       and that all refernces to protons will become neutrons if that is the incoming nucleon
    // NOTE: this class works for all neutrino flavours

    class pv_pv : public ME2_Base
    {
    private:
    public:
        pv_pv(const External_ME_Args &args, const incomingnucleon::code &nucleon);

        double operator()(const ATOOLS::Vec4D_Vector &mom);
        double m_alpha, m_alphas, m_s, m_massZ;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
        incomingnucleon::code m_nucleon;
    };

    pv_pv::pv_pv(const External_ME_Args &args, const incomingnucleon::code &nucleon)
        : ME2_Base(args), m_nucleon(nucleon)
    {
        msg_Out() << "pv_pv::pv_pv(): Constructor called" << std::endl;
        m_sintt = 2; // T-channel only , maybe look at single channel restriction later
        // m_sintt = 0;  // No s/t/u channel restriction
        m_oew = 0;  // EW order zero (no loops)
        m_oqcd = 0; // QCD order zero
        msg_Out() << "pv_pv::pv_pv(): Getting constants" << std::endl;
        // m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");
        // m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");
        m_massZ = Flavour(kf_Z).Mass();

        msg_Out() << "pv_pv::pv_pv(): Creating FormFactor_EMnucleon for Z" << std::endl;
        p_formfactor = std::make_unique<FormFactor_EMnucleon>(incomingboson::Z, m_nucleon);

        msg_Out() << "pv_pv::pv_pv(): Constructor complete for nucleon " << m_nucleon << std::endl;
    }

    double pv_pv::operator()(const ATOOLS::Vec4D_Vector &momenta)
    {
        // indices: # p[0] v[1] -> p[2] v[3]
        // k: lepton, p: hadron
        const auto &ki = momenta[1];
        const auto &pi = momenta[0];
        const auto &kf = momenta[3];
        const auto &pf = momenta[2];

        // masses squared (on shell)
        const double Mp2 = pi.Abs2();

        // Mandelstraam variables
        double s = (ki + pi).Abs2();
        double t = (ki - kf).Abs2(); // t < 0 for spacelike
        double u = (ki - pf).Abs2();

        double Q2 = -t; // Q2 > 0 for spacelike photon exchange

        // Add a small cutoff to avoid IR divergence at Q2=0 (forward scattering)
        // cutoff order approx mass of electron, we dont expect slow moving electrons
        const double Q2_min = 1e-3;
        if (Q2 < Q2_min)
        {
            // msg_Out() << "pv_pv::operator(): Q2 < Q2_min, returning 0.0" << std::endl;
            return 0.0;
        }

        NucleonFormFactors ff = p_formfactor->GetFormFactors(Q2);
        double F1 = ff.F1;
        double F2 = ff.F2;
        double FA_Z = ff.FA;
        double FP_Z = ff.FP;
        double F12 = F1 + F2;

        double A = (s - Mp2) / 2.0;  // (ki . pi)
        double B = -(u - Mp2) / 2.0; // (ki . pf)
        double C = -(t) / 2.0;       // (ki . kf)
        double D = Mp2 - (t) / 2.0;  // (pi . pf)

        // Calculate spin averaged matrix element squared
        double term1coeff = 4 * F12 * F12;
        double term1 = term1coeff * (2 * A * A + 2 * B * B + C * t + D * t);
        double term2coeff = 2 * ((F2 * F2 * (Mp2 + 4 * D)) / (8 * Mp2) - 2 * F2 * F12);
        double term2 = term2coeff * (2 * (A * A + B * B + 2 * A * B) + (C + D + 2 * B) * t * 0.5);

        double L_munu_H_munu = -(term1 + term2);
        double pre_factor = 32 * M_PI * M_PI * sqr((*aqed)(m_s)) / (Q2 * Q2);

        double M_squared = pre_factor * L_munu_H_munu;

        return M_squared;
    }
}

DECLARE_TREEME2_GETTER(EXTRAXS::pv_pv, "pv_pv")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::pv_pv>::
operator()(const External_ME_Args &args) const
{
    const Flavour_Vector &fl = args.Flavours();
    incomingnucleon::code nucleon = incomingnucleon::off;

    // check if elastic scattering
    if (fl.size() != 4)
        return NULL;
    if (fl[0] != fl[2])
        return NULL;
    if (fl[1] != fl[3])
        return NULL;

    // Sherpa orders: fl[0]=hadron_in, fl[1]=lepton_in, fl[2]=hadron_out, fl[3]=lepton_out
    if (fl[0] == Flavour(kf_p_plus))
    {
        nucleon = incomingnucleon::proton;
    }
    else if (fl[0] == Flavour(kf_n))
    {
        nucleon = incomingnucleon::neutron;
    }
    else
    {
        return NULL; // Not a nucleon
    }

    // Check for neutrino (NC interaction)
    if (fl[1].IsNeutrino())
    {
        return new pv_pv(args, nucleon);
    }

    return NULL;
}
