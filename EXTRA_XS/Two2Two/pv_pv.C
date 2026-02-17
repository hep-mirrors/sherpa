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
        pv_pv(const External_ME_Args &args, const incomingnucleon::code &nucleon, 
            const Flavour &neutrino_fl, const Flavour &nucleon_fl);

        double operator()(const ATOOLS::Vec4D_Vector &mom);
        double m_massZ;
        double m_sin2_theta_w;
        double m_g;
        double m_cos_theta_w;
        double m_gz;
        Flavour m_fl_neutrino, m_fl_nucleon;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
        incomingnucleon::code m_nucleon;
    };

    pv_pv::pv_pv(const External_ME_Args &args, const incomingnucleon::code &nucleon, 
        const Flavour &neutrino_fl, const Flavour &nucleon_fl)
        : ME2_Base(args), m_nucleon(nucleon), m_fl_neutrino(neutrino_fl), m_fl_nucleon(nucleon_fl)
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
        m_sin2_theta_w = 0.23126;
        m_g = 0.653;
        m_cos_theta_w = sqrt(1.0 - m_sin2_theta_w);
        m_gz = m_g / m_cos_theta_w;

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

        msg_Out() << "pe_nv::operator(): ki = " << ki << std::endl;
        msg_Out() << "pe_nv::operator(): pi = " << pi << std::endl;
        msg_Out() << "pe_nv::operator(): kf = " << kf << std::endl;
        msg_Out() << "pe_nv::operator(): pf = " << pf << std::endl;

        // Calculate scattering angle theta between initial and final lepton
        Vec3D ki_vec(ki);
        Vec3D kf_vec(kf);
        double cos_theta = (ki_vec * kf_vec) / (ki_vec.Abs() * kf_vec.Abs());
        double theta = acos(cos_theta);
        double phi = atan2(kf[2], kf[1]);

        msg_Out() << "pv_pv::operator(): theta = " << theta << std::endl;
        msg_Out() << "pv_pv::operator(): phi = " << phi << std::endl;

        // masses squared (on shell)
        const double m = m_fl_nucleon.Mass();
        const double m2 = m * m;
        msg_Out() << "pv_pv::operator():  m2 = " << m2 << std::endl;

        // Mandelstraam variables
        double s = (ki + pi).Abs2();
        double t = (ki - kf).Abs2(); // t < 0 for spacelike
        double u = (ki - pf).Abs2();

        double Q2 = -t; // Q2 > 0 for spacelike photon exchange
        // msg_Out() << "pe_pe::operator(): Q2 = " << Q2 << std::endl;
        // msg_Out() << "pe_pe::operator(): s = " << s << ", t = " << t << ", u = " << u << std::endl;

        // Add a small cutoff to avoid IR divergence at Q2=0 (forward scattering)
        // cutoff order less than mass of electron
        const double Q2_min = 1e-5;
        if (Q2 < Q2_min)
        {
            // msg_Out() << "pv_pv::operator(): Q2 < Q2_min, returning 0.0" << std::endl;
            return 0.0;
        }

        NucleonFormFactors ff = p_formfactor->GetFormFactors(Q2);
        const double F1 = ff.F1;
        const double F2 = ff.F2;
        const double FA_Z = ff.FA;
        const double FP_Z = ff.FP;

        // Get coupling constants and masses
        double gA = 0.5; // Axial coupling (Z boson)
        double gV = -0.5; // Vector coupling (Z boson)
        double FA = FA_Z;

        // Precalculate 
        const double MZ2 = m_massZ * m_massZ;
        const double s2 = s*s;
        const double u2 = u*u;
        const double sum1 = gA *gA + gV *gV;
        const double sum2 = s-t+u;
        const double m4 = m2 * m2;
        const double gz4 = m_gz * m_gz * m_gz * m_gz;

        // calculate
        double term1 = 2.0 * F1 * F1 * (sum1) *
                        (2.0 * m4 - 2.0 * m2 * (sum2) + s2 + u2);

        double term2 = 4.0 * F1 * (2.0 * m2 - s - u) * 
                       (F2 * t * (sum1) + 2.0 * FA * gA * gV * (u-s));
        
        double term3 = F2 * F2 * t * (sum1) * 
                       (m2 - (sum2) + (1/m2) * s * u);
        
        double term4 = 8.0 * F2 * FA * gA * gV * t * (u-s);
        
        double term5 = 2.0 * FA * FA * (sum1) * 
                       (2.0 * m4 - 2.0 * m2 * (s+t+u) + s2 + u2);
        
        double coeff = gz4 / ((MZ2 + Q2) * (MZ2 + Q2));
        double M_squared = coeff * 
                          (term1 + term2 + term3 + term4 + term5);
        
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
    // if (fl.size() != 4)
    //     return NULL;
    msg_Out() << "pv_pv Getter: Flavour vector size = " << fl.size() << std::endl;
    for (size_t i = 0; i < fl.size(); ++i)
    {
        msg_Out() << "pv_pv Getter: fl[" << i << "] = " << fl[i] << std::endl;
    }
    // Sherpa orders: fl[0]=hadron_in, fl[1]=lepton_in, fl[2]=hadron_out, fl[3]=lepton_out ??
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
        return new pv_pv(args, nucleon, fl[1], fl[0]);
    }

    return NULL;
}
