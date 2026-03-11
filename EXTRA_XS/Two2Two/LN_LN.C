#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"
#include <memory>
#include <stdexcept>

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;
using namespace PHASIC;
using namespace std;

// ============================================================
//  Assumed momentum ordering: 
//  (use Sort_Flavours=0 in runcard under PROCESSES to allow this)
//
//    momenta[0] = lepton_in   (e-, mu-, ve, ...)
//    momenta[1] = nucleon_in  (p, n)
//    momenta[2] = lepton_out  (e-, mu-, ve, ...)
//    momenta[3] = nucleon_out (p,n)
// ============================================================

namespace EXTRAXS
{
    // ==========================================================
    //  SHARED BASE CLASS
    //  Holds: incoming nucleon/lepton flavours, kinematic helper, Q2 cut
    //  TODO: rewrite form factors to remove incomingnucleon code and depend only on flavour
    // ==========================================================
    class LN_LN_Base : public ME2_Base
    {
    protected:
        incomingnucleon::code m_nucleon;
        Flavour m_fl_lepton;
        Flavour m_fl_nucleon;

        static constexpr double s_Q2min = 1e-5;

        struct Kinematics {
            ATOOLS::Vec4D ki, pi, kf, pf;   // individual momenta for checking, will remove in final versions
            double s, t, u, Q2;
        };
 
        Kinematics ComputeKinematics(const ATOOLS::Vec4D_Vector &mom) const
        {
            Kinematics k;
            k.ki = mom[0];  k.pi = mom[1];
            k.kf = mom[2];  k.pf = mom[3];
            k.s  = (k.ki + k.pi).Abs2();
            k.t  = (k.ki - k.kf).Abs2();   // t < 0 (spacelike)
            k.u  = (k.ki - k.pf).Abs2();
            k.Q2 = -k.t;                   // Q2 > 0
            return k;
        }

    public:
        LN_LN_Base(const External_ME_Args &args,
                   const incomingnucleon::code &nucleon,
                   const Flavour &lepton_fl,
                   const Flavour &nucleon_fl)
            : ME2_Base(args),
              m_nucleon(nucleon),
              m_fl_lepton(lepton_fl),
              m_fl_nucleon(nucleon_fl)
        {
            m_sintt = 2;   // t-channel topology only
            m_oew   = 0;   // tree-level EW
            m_oqcd  = 0;   // no QCD
        }

        virtual double operator()(const ATOOLS::Vec4D_Vector &mom) = 0;
        virtual ~LN_LN_Base() = default;
    };


    // ==========================================================
    //  NC CHARGED LEPTON:  l N -> l N   (photon [+ Z])
    // ==========================================================
    class LN_LN_nc_charged : public LN_LN_Base
    {
    private:
        std::unique_ptr<FormFactor_EMnucleon> p_ff_photon;
        std::unique_ptr<FormFactor_EMnucleon> p_ff_Z;      // only if m_include_Z
        double m_alpha, m_e2, m_massZ;
        bool   m_include_Z;

    public:
        LN_LN_nc_charged(const External_ME_Args &args,
                    const incomingnucleon::code &nucleon,
                    const Flavour &lepton_fl,
                    const Flavour &nucleon_fl)
            : LN_LN_Base(args, nucleon, lepton_fl, nucleon_fl)
        {
            msg_Out() << "\n\n"
                      << "----------------------------------------------\n"
                      << METHOD << ": e N -> e N (NC charged lepton)" << "\n";

            m_alpha     = 1.0 / 137.0;
            m_e2        = 4.0 * M_PI * m_alpha;
            m_massZ     = Flavour(kf_Z).Mass();

            ATOOLS::Settings &s = ATOOLS::Settings::GetMainSettings();
            m_include_Z = s["Form_Factor"]["Include_Z_Interference"].SetDefault(true).Get<bool>();

            msg_Out() << "  -> initializing photon form factors for incoming nucleon: " << m_nucleon << "\n";
            p_ff_photon = std::make_unique<FormFactor_EMnucleon>(incomingboson::photon, m_nucleon);
            if (m_include_Z)
                msg_Out() << "  -> initializing Z form factors for incoming nucleon: " << m_nucleon << "\n";
                p_ff_Z = std::make_unique<FormFactor_EMnucleon>(incomingboson::Z, m_nucleon);
        }

        double operator()(const ATOOLS::Vec4D_Vector &mom) override
        {
            Kinematics k = ComputeKinematics(mom);
            if (k.Q2 < s_Q2min) return 0.0;
            // Mandelstraam variables
            const double s = k.s;
            const double t = k.t; // t < 0 virtual boson exchange
            const double u = k.u; // u very small, so dont use in calculation to avoid rounding errors.
            const double Q2 = k.Q2;

            NucleonFormFactors ff = p_ff_photon->GetFormFactors(k.Q2);
            const double F1 = ff.F1;
            const double F2 = ff.F2;
            // msg_Out() << "pe_pe::operator(): F1 = " << F1 << ", F2 = " << F2 << std::endl;

            const double Me2 = Flavour(kf_e).Mass()*Flavour(kf_e).Mass();
            const double Mp2 = k.pi.Abs2();  
            
            // Pre-compute common combinations
            const double s2 = s * s;
            const double u2 = u * u;
            const double t2 = t * t;
            const double Me4 = Me2 * Me2;
            const double Me6 = Me4 * Me2;
            const double Mp4 = Mp2 * Mp2;
            const double F1_2 = F1 * F1;
            const double F2_2 = F2 * F2;
            const double twoMe2 = 2 * Me2;
            const double twoMp2 = 2 * Mp2;
            const double twoMe2_plus_t = twoMe2 + t;
            const double twoMe2_plus_twoMp2_minus_s_minus_u = twoMe2 + twoMp2 - s - u;

            double calc =
                2 * F1_2 * Mp2 * (2 * Me4 + twoMe2 * (twoMp2 - s + t - u) + 2 * Mp4 - twoMp2 * (s - t + u) + s2 + u2) + 4 * F1 * F2 * Mp2 * twoMe2_plus_t * twoMe2_plus_twoMp2_minus_s_minus_u + F2_2 * (4 * Me6 + Me4 * (8 * Mp2 - 4 * s + t - 4 * u) + Me2 * (4 * Mp4 + Mp2 * (-4 * s + 6 * t - 4 * u) + s2 - s * t + 2 * s * u - t2 - t * u + u2) + t * (Mp4 - Mp2 * (s - t + u) + s * u));

            const double Q4 = Q2 * Q2;
            const double pre_factor = m_e2 * m_e2 / (Mp2 * Q4);
            double M_squared = pre_factor * calc;

            if (m_include_Z) {
                NucleonFormFactors ff_Z = p_ff_Z->GetFormFactors(k.Q2);
                // TODO: paste Z-interference calculation from pe_pe.C here
                M_squared += 0.0; // placeholder for Z contribution and interference term
            }

            return M_squared;
        }
    };


    // ==========================================================
    //  NC NEUTRINO:  v N -> v N   (Z only)
    // ==========================================================
    class LN_LN_nc_neutrino : public LN_LN_Base
    {
    private:
        std::unique_ptr<FormFactor_EMnucleon> p_ff_Z;
        double m_massZ, m_g, m_sin2_theta_w, m_cos_theta_w, m_gz;

    public:
        LN_LN_nc_neutrino(const External_ME_Args &args,
                          const incomingnucleon::code &nucleon,
                          const Flavour &neutrino_fl,
                          const Flavour &nucleon_fl)
            : LN_LN_Base(args, nucleon, neutrino_fl, nucleon_fl)
        {
            msg_Out() << "\n\n"
                      << "----------------------------------------------\n"
                      << METHOD << ": v N -> v N (NC neutrino)" << "\n";

            m_massZ        = Flavour(kf_Z).Mass();
            m_sin2_theta_w = 0.23126;
            m_g            = 0.653;
            m_cos_theta_w  = sqrt(1.0 - m_sin2_theta_w);
            m_gz           = m_g / m_cos_theta_w;

            msg_Out() << "  -> initializing Z form factors for incoming nucleon: " << m_nucleon << "\n";
            p_ff_Z = std::make_unique<FormFactor_EMnucleon>(incomingboson::Z, m_nucleon);
        }

        double operator()(const ATOOLS::Vec4D_Vector &mom) override
        {
            Kinematics k = ComputeKinematics(mom);
            if (k.Q2 < s_Q2min) return 0.0;

            NucleonFormFactors ff = p_ff_Z->GetFormFactors(k.Q2);
            // TODO: paste calculation from pv_pv.C here

            return 0.0;
        }
    };

    // ==========================================================
    //  CC:  l N -> v N'  &  v N' -> l N    (W)
    // ==========================================================
    class LN_LN_cc : public LN_LN_Base
    {
    private:
        std::unique_ptr<FormFactor_EMnucleon> p_ff_W;
        double m_massW, m_g, m_gA;

    public:
        LN_LN_cc(const External_ME_Args &args,
                    const incomingnucleon::code &nucleon,
                    const Flavour &lepton_fl,
                    const Flavour &nucleon_fl)
            : LN_LN_Base(args, nucleon, lepton_fl, nucleon_fl)
        {
            msg_Out() << "\n\n"
                      << "----------------------------------------------\n"
                      << METHOD << ": l N -> l' N' (CC)" << "\n";

            m_massW = Flavour(kf_Wplus).Mass();
            m_g     = 0.653;
            m_gA    = 1.27;

            msg_Out() << "  -> initializing W form factors for incoming nucleon: " << m_nucleon << "\n";
            p_ff_W = std::make_unique<FormFactor_EMnucleon>(incomingboson::W, m_nucleon);
        }

        double operator()(const ATOOLS::Vec4D_Vector &mom) override
        {
            Kinematics k = ComputeKinematics(mom);
            if (k.Q2 < s_Q2min) return 0.0;

            NucleonFormFactors ff = p_ff_W->GetFormFactors(k.Q2);
            // masses (on shell)
            // Mass variables (not squared)
            const double Me = Flavour(kf_e).Mass();
            const double Mp = Flavour(kf_p_plus).Mass();
            const double Mn = Flavour(kf_n).Mass();
            const double Me2 = Me * Me;
            const double Mp2 = Mp * Mp;
            const double Mn2 = Mn * Mn;
            msg_Out() << "pe_nv::operator(): Me2 = " << Me2 << ", Mp2 = " << Mp2 << ", Mn2 = " << Mn2 << std::endl;
            // msg_Out() << " inv mass from mom: m_ki^2" << ki.Abs2() << "m_pi ^ 2 " << pi.Abs2() << "m_kf^2" << kf.Abs2() << "m_pf^2" << pf.Abs2() << std::endl;
            // Mandelstraam variables
            const double s = k.s;
            const double t = k.t; // t < 0 virtual boson exchange
            const double u = k.u; // u very small, so dont use in calculation to avoid rounding errors.

            const double Q2 = k.Q2; // Q2 > 0 for virtual boson exchange
            msg_Out() << "pe_nv::operator(): Q2 = " << Q2 << std::endl;
            msg_Out() << "pe_nv::operator(): s = " << s << ", t = " << t << ", u = " << u << std::endl;

            // Get form factors (constants for given Q2)
            NucleonFormFactors ff = p_ff_W->GetFormFactors(Q2);
            const double F1 = ff.F1;
            const double F2 = ff.F2;
            const double FA = ff.FA;
            const double FP = ff.FP;
            msg_Out() << "pe_nv::operator(): F1 = " << F1 << ", F2 = " << F2 << ", FA = " << FA << ", FP" << FP << std::endl;

            // precalcs
            const double MW2 = m_massW * m_massW;
            const double gW4 = (m_g * m_g * m_g * m_g) / 4;
            const double sum_mass = Mn + Mp;
            const double sum_mass2 = (Mn + Mp) * (Mn + Mp);
            const double diff_mass = Mn - Mp;
            const double diff_mass2 = (Mn - Mp) * (Mn - Mp);
            const double t2 = t * t;
            const double s2 = s * s;
            const double st = s * t;

            // Calc
            double calc = Mp2 * (
                F1 * F1 * (
                    Me2 * (diff_mass2 - 2.0 * s - t) + 2.0 * (Mn2 - s) * (Mp2 - s) - t * diff_mass2 + 2.0 * st + t2
                ) 
                + 2.0 * F1 * FA * (
                    Me2 * (Mn2 - Mp2 + t) + t * (Mn2 + Mp2 - 2.0 * s - t)
                ) 
                + FA * FA * (
                    Me2 * (sum_mass2 - 2.0 * s - t) + 2.0 * (Mn2 - s) * (Mp2 - s) - t * sum_mass2 + 2.0 * st + t2
                )
            ) 
            - F2 * Mp * (
                    F1 * (Me2 * Me2 * (-Mn) + Me2 * (-(diff_mass * (Mp2 - s)) - Mp * t) - t * sum_mass * (diff_mass2 - t)) 
                    + FA * sum_mass * (Me2 * (Mn2 - Mp2 + t) + t * (Mn2 + Mp2 - 2.0 * s - t))
                ) 
            + (1.0 / 8.0) * F2 * F2 * (
                -(Me2 * Me2 * ((3.0 * Mn - Mp) * sum_mass + t)) + Me2 * (-2.0 * diff_mass * sum_mass * (Mn2 + Mp2 - 2.0 * s) + t * (Mn - 3.0 * Mp) * sum_mass + 4.0 * st + t2) - 2.0 * t * (Mn2 * Mn2 - 2.0 * s * (Mn2 + Mp2) + Mp2 * Mp2 + 2.0 * s2) + 2.0 * t2 * (sum_mass2 - 2.0 * s)
            ) 
            + FA * FP * Me2 * Mp * (
                Me2 * Mn + Mn2 * Mn - Mn2 * Mp - Mn * (s + t) + Mp * s
            ) 
            + (1.0 / 8.0) * FP * FP * Me2 * (Me2 - t) * (diff_mass2 - t);

            // Coefficient
            const double prefactor = gW4 / (4.0 * Mp2 * (MW2 + Q2) * (MW2 + Q2));

            double M_squared = prefactor * calc;

            return M_squared;
        }
    };

} // namespace EXTRAXS


// ============================================================
//  SINGLE GETTER  — dispatches all three sub-processes
//
//  Expected ordering:
//    fl[0] = lepton/neutrino in   fl[1] = nucleon in
//    fl[2] = lepton/neutrino out  fl[3] = nucleon out
// ============================================================
DECLARE_TREEME2_GETTER(EXTRAXS::LN_LN_Base, "LN-LN")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::LN_LN_Base>::
operator()(const External_ME_Args &args) const
{
    const Flavour_Vector &fl = args.Flavours();
    if (fl.size() != 4) return NULL;

    // Require IS: lepton + nucleon
    if (!fl[0].IsLepton()  || !fl[1].IsNucleon()) return NULL;
    // Require FS: lepton + nucleon
    if (!fl[2].IsLepton()  || !fl[3].IsNucleon()) return NULL;

    msg_Out() << "\n----------------------------------------------\n"
              << METHOD << ": "
              << fl[0] << " " << fl[1] << " --> " << fl[2] << " " << fl[3] << "\n";

    // Identify incoming nucleon
    incomingnucleon::code nucleon = incomingnucleon::off;
    if      (fl[1] == Flavour(kf_p_plus)) nucleon = incomingnucleon::proton;
    else if (fl[1] == Flavour(kf_n))      nucleon = incomingnucleon::neutron;
    else    return NULL;

    // ---- CC: charged lepton -> neutrino (or vice-versa) ----
    // Lepton identity changes between IS and FS
    if (fl[0].IsChargedLepton() && fl[2].IsNeutrino()) {
        msg_Out() << "  -> dispatching to LN_LN_cc \n";
        return new EXTRAXS::LN_LN_cc(args, nucleon, fl[0], fl[1]);
    }
    if (fl[0].IsNeutrino() && fl[2].IsChargedLepton()) {
        // Reverse CC: nu N -> l N'  — nucleon identity also reverses
        msg_Out() << "  -> dispatching to LN_LN_cc \n";
        return new EXTRAXS::LN_LN_cc(args, nucleon, fl[0], fl[1]);
    }

    // ---- NC: lepton identity preserved (fl[0]==fl[2], fl[1]==fl[3]) ----
    if (fl[0] != fl[2] || fl[1] != fl[3]) return NULL;

    if (fl[0].IsChargedLepton()) {
        msg_Out() << "  -> dispatching to LN_LN_nc_charged \n";
        return new EXTRAXS::LN_LN_nc_charged(args, nucleon, fl[0], fl[1]);
    }
    if (fl[0].IsNeutrino()) {
        msg_Out() << "  -> dispatching to LN_LN_nc_neutrino \n";
        return new EXTRAXS::LN_LN_nc_neutrino(args, nucleon, fl[0], fl[1]);
    }

    return NULL;
}
