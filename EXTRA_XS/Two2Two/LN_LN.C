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

        static constexpr double s_Q2min = 1e-5;

        struct Kinematics {
            ATOOLS::Vec4D ki, pi, kf, pf;   // individual momenta for checking, will remove in final versions
            double s, t, u, Q2;
            double theta, phi;              // for checking momentum input
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
            Vec3D ki_vec(k.ki);
            Vec3D kf_vec(k.kf);
            double cos_theta = (ki_vec * kf_vec) / (ki_vec.Abs() * kf_vec.Abs());
            k.theta = acos(cos_theta);
            k.phi   = atan2(k.kf[2], k.kf[1]);
            return k;
        }

    public:
        LN_LN_Base(const External_ME_Args &args,
                   const incomingnucleon::code &nucleon)
            : ME2_Base(args),
              m_nucleon(nucleon)
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
        double m_alpha, m_Qf, m_e2, m_massZ2, m_gZ4, m_g, m_sin2thetaW, m_costhetaW, m_cV, m_cA;
        double m_mass2L, m_mass2N, m_mass4L, m_mass4N;
        bool   m_include_Z;

    public:
        LN_LN_nc_charged(const External_ME_Args &args,
                    const incomingnucleon::code &nucleon)
            : LN_LN_Base(args, nucleon)
        {
            const Flavour_Vector &fl = args.Flavours();
            msg_Out() << "\n\n"
                      << "----------------------------------------------\n"
                      << METHOD << ": e N -> e N (NC charged lepton)" << "\n"
                      << "check flavours " << fl[0] << " " << fl[1] << " -> " << fl[2] << " " << fl[3] << "\n";


            m_alpha     = 1.0 / 137.0;
            m_e2        = 4.0 * M_PI * m_alpha;
            m_Qf        = -1;
            m_mass2L    = fl[0].Mass() * fl[0].Mass();
            m_mass2N    = fl[1].Mass() * fl[1].Mass();
            m_mass4L    = m_mass2L*m_mass2L;
            m_mass4N    = m_mass2N * m_mass2N;

            ATOOLS::Settings &s = ATOOLS::Settings::GetMainSettings();
            m_include_Z = s["Form_Factor"]["Include_Z_Interference"].SetDefault(true).Get<bool>();

            msg_Out() << "  -> initializing photon form factors for incoming nucleon: " << m_nucleon << "\n";
            p_ff_photon = std::make_unique<FormFactor_EMnucleon>(incomingboson::photon, m_nucleon);
            if (m_include_Z)
                msg_Out() << "  -> initializing Z form factors for incoming nucleon: " << m_nucleon << "\n";
                p_ff_Z = std::make_unique<FormFactor_EMnucleon>(incomingboson::Z, m_nucleon);
                m_massZ2     = Flavour(kf_Z).Mass()*Flavour(kf_Z).Mass();
                m_sin2thetaW = 0.23126;
                m_costhetaW = sqrt(1.0 - m_sin2thetaW);
                m_g = 0.653;
                m_gZ4 = pow(m_g / m_costhetaW, 4);
                m_cA = -0.5;
                m_cV = -0.5 + 2*m_sin2thetaW;
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

            // Pre-compute common combinations
            const double me2 = m_mass2L;
            const double mp2 = m_mass2N;
            const double me4 = m_mass4L;
            const double mp4 = m_mass4N;

            const double delta = me2 + mp2 - s;
            const double delta2 = delta * delta;
            const double t2 = t * t;
            const double st = s * t;
            const double s_plus_t = s + t;

            const double calc = 2.0 * F1*F1 * mp2 * (2.0 * delta2 + 2.0 * st + t2) 
                + 4.0 * F1 * F2 * mp2 * t * (2.0 * me2 + t) 
                + F2*F2 * t * (-me4 + me2 * (2.0 * (mp2 + s) + t) - mp4 + 2.0 * mp2 * s_plus_t - s * s_plus_t);

            const double pre_factor = m_e2 * m_e2 / (m_mass2N * Q2*Q2);
            double M_squared = pre_factor * calc;

            if (m_include_Z) {
                NucleonFormFactors ff_Z = p_ff_Z->GetFormFactors(k.Q2);
                const double F1_Z = ff.F1;
                const double F2_Z = ff.F2;
                const double FA_Z = ff.FA;
                const double FP_Z = ff.FP;

                const double mp2_minus_s = mp2 - s;
                const double mp2_minus_s2 = mp2_minus_s * mp2_minus_s;
                const double s2 = s*s;

                // cA^2 block
                const double CA_term1 =
                    2.0 * F1_Z*F1_Z * mp2 *
                    (2.0 * me4 - 4.0 * me2 * (mp2 + s + t) + 2.0 * mp2_minus_s2 + 2.0 * st + t2);

                const double CA_term2 =
                    4.0 * F1_Z * F2_Z * mp2 * t * (t - 4.0 * me2);

                const double CA_term3 =
                    t * (FP_Z*FP_Z * me2 * t - F2_Z*F2_Z * (me4 + me2 * (6.0 * mp2 - 2.0 * s) + mp4 - 2.0 * mp2 * s_plus_t + s * s_plus_t));

                const double CA_term4 =
                    2.0 * FA_Z*FA_Z * mp2 *
                    (2.0 * me4 + 4.0 * me2 * (5.0 * mp2 - s - t) + 2.0 * mp4 - 4.0 * mp2 * s_plus_t + 2.0 * s2 + 2.0 * st + t2);

                const double CA_term5 =
                    8.0 * FA_Z * FP_Z * me2 * mp2 * t;

                const double CA_block = m_cA*m_cA * (CA_term1 + CA_term2 + CA_term3 + CA_term4 + CA_term5);

                // mixed cA cV block
                const double CAcV_block =
                    8.0 * m_cA*m_cV * FA_Z * mp2 * t * (F1_Z + F2_Z) * (2.0 * delta - t);

                // cV^2 block
                const double CV_term1 =
                    2.0 * F1_Z*F1_Z * mp2 * (2.0 * delta2 + 2.0 * st + t2);

                const double CV_term2 =
                    4.0 * F1_Z * F2_Z * mp2 * t * (2.0 * me2 + t);

                const double CV_term3 =
                    F2_Z*F2_Z * t *
                    (-me4 + me2 * (2.0 * (mp2 + s) + t) - mp4 + 2.0 * mp2 * s_plus_t - s * s_plus_t);

                const double CV_term4 =
                    2.0 * FA_Z*FA_Z * mp2 *
                    (2.0 * me4 - 4.0 * me2 * (mp2 + s) + 2.0 * mp4 - 4.0 * mp2 * s_plus_t + 2.0 * s2 + 2.0 * st + t2);

                const double CV_block = m_cV*m_cV * (CV_term1 + CV_term2 + CV_term3 + CV_term4);

                // final Z calc
                const double calc_Z = CA_block + CAcV_block + CV_block;

                const double prefactor_Z = m_gZ4 / (4*mp2*(m_massZ2 + Q2)*(m_massZ2 + Q2));
                double M2_Z = prefactor_Z*calc_Z;
                double M2 = M_squared;

                M_squared += 2*sqrt(M2)*sqrt(M2_Z) + M2_Z ; // placeholder for Z contribution and interference term
            }

            return M_squared;
        }
    };


    // ==========================================================
    //  NC NEUTRINO:  v N -> v N   (Z only)
    //      TODO: check couplings, otherwise correct.
    // ==========================================================
    class LN_LN_nc_neutrino : public LN_LN_Base
    {
    private:
        std::unique_ptr<FormFactor_EMnucleon> p_ff_Z;
        double m_massZ2, m_g, m_gA, m_gV, m_sin2_theta_w, m_cos_theta_w, m_gZ4;
        double m_mass2N, m_mass4N;
        double m_helicityavg;

    public:
        LN_LN_nc_neutrino(const External_ME_Args &args,
                          const incomingnucleon::code &nucleon)
            : LN_LN_Base(args, nucleon)
        {
            const Flavour_Vector &fl = args.Flavours();
            msg_Out() << "\n\n"
                      << "----------------------------------------------\n"
                      << METHOD << ": v N -> v N (NC neutrino)" << "\n"
                      << "check flavours " << fl[0] << " " << fl[1] << " -> " << fl[2] << " " << fl[3] << "\n";

            m_massZ2       = Flavour(kf_Z).Mass()*Flavour(kf_Z).Mass();
            m_sin2_theta_w = 0.23126;
            m_g            = 0.653;
            m_gA           = 0.5;
            m_gV           = 0.5;
            m_cos_theta_w  = sqrt(1.0 - m_sin2_theta_w);
            m_gZ4          = pow(m_g / m_cos_theta_w,4);
            m_mass2N       = fl[1].Mass() * fl[1].Mass();
            m_mass4N       = m_mass2N * m_mass2N;
            m_helicityavg  = 1;

            msg_Out() << "  -> initializing Z form factors for incoming nucleon: " << m_nucleon << "\n";
            p_ff_Z = std::make_unique<FormFactor_EMnucleon>(incomingboson::Z, m_nucleon);
        }

        double operator()(const ATOOLS::Vec4D_Vector &mom) override
        {
            Kinematics k = ComputeKinematics(mom);
            if (k.Q2 < s_Q2min) return 0.0;

            const double s = k.s;
            const double t = k.t; // t < 0 virtual boson exchange
            const double u = k.u; // u very small, so dont use in calculation to avoid rounding errors.

            const double Q2 = k.Q2; // Q2 > 0 for virtual boson exchange
            msg_Debugging() << METHOD << ": theta = " << k.theta << ", phi = " << k.phi << std::endl;
            msg_Debugging() << METHOD << ": Q2 = " << Q2 << std::endl;
            msg_Debugging() << METHOD << ": s = " << s << ", t = " << t << ", u = " << u << std::endl;
            msg_Debugging() << METHOD << ": ki = " << k.ki << ", pi = " << k.pi << ", kf = " << k.kf << ", pf = " << k.pf <<"\n";
            msg_Debugging() << METHOD << ": MassNucleon^2 = " << m_mass2N << "\n" 
                            << METHOD << ": MassZ^2       = " << m_massZ2 << "\n";

            NucleonFormFactors ff = p_ff_Z->GetFormFactors(k.Q2);
            const double F1 = ff.F1;
            const double F2 = ff.F2;
            const double FA = ff.FA;
            msg_Debugging() << METHOD << ": F1 = " << F1 << ", F2 = " << F2 << ", FA = " << FA << ", FP not needed."  
                            << METHOD << ": gZ^4 = " <<  m_gZ4 << std::endl;

            // precalcs
            const double m2 = m_mass2N;
            const double m4 = m_mass4N;
            const double s2 = s * s;
            const double u2 = u * u;
            const double su = s * u;
            const double u_minus_s = u - s;
            const double kin_stu = s - t + u;
            const double kin_plus = s + t + u;

            // repeated kinematic blocks
            const double K1 = 2.0 * m4 - 2.0 * m2 * kin_stu + s2 + u2;
            const double K2 = 2.0 * m4 - 2.0 * m2 * kin_plus + s2 + u2;
            const double K3 = m4 - m2 * kin_stu + su;

            const double gA2_plus_gV2 = m_gA*m_gA + m_gV*m_gV;

            const double A1 = -2.0 * m2 + s + u;
            const double A2 = F2 * t * gA2_plus_gV2 + 2.0 * FA * m_gA*m_gV * u_minus_s;

            // calc
            const double calc =
                2.0 * F1*F1 * m2 * gA2_plus_gV2 * K1
                - 4.0 * F1 * m2 * A1 * A2
                + F2*F2 * t * gA2_plus_gV2 * K3
                + 8.0 * F2 * FA * m_gA*m_gV * m2 * t * u_minus_s
                + 2.0 * FA*FA * m2 * gA2_plus_gV2 * K2;

            const double prefactor = (4 * m_gZ4) / (m_helicityavg * m2 * (m_massZ2 + Q2)*(m_massZ2 + Q2));

            double M2 = prefactor * calc;
            msg_Debugging() << METHOD << ": M2 = " << M2 << std::endl;

            return M2;
        }
    };

    // ==========================================================
    //  CC:  v N -> l N'  &  l N' -> v N    (W)
    // ==========================================================
    class LN_LN_cc : public LN_LN_Base
    {
    private:
        std::unique_ptr<FormFactor_EMnucleon> p_ff_W;
        double m_massWsqr, m_g, m_gW4, m_gA;
        double m_massL, m_massp, m_massn;
        double m_helicityavg, m_massavg;

    public:
        LN_LN_cc(const External_ME_Args &args,
                    const incomingnucleon::code &nucleon, bool incomingneutrino)
            : LN_LN_Base(args, nucleon)
        {
            const Flavour_Vector &fl = args.Flavours();
            msg_Out() << "\n\n"
                      << "----------------------------------------------\n"
                      << METHOD << ": l N -> l' N' (CC)" << "\n"
                      << "check flavours " << fl[0] << " " << fl[1] << " -> " << fl[2] << " " << fl[3] << "\n";

            m_massWsqr = Flavour(kf_Wplus).Mass() * Flavour(kf_Wplus).Mass();
            m_g     = 0.653;
            m_gW4    = (m_g*m_g*m_g*m_g)/4;
            m_gA    = 1.27;
            if (incomingneutrino == true){
                m_helicityavg = 2; 
            } else if (incomingneutrino == false){
                m_helicityavg = 4;
            }
            else {THROW(missing_input,"Missing incoming neutrino/lepton beam information.")}
            const double incoming_lepton_mass = fl[0].Mass();
            m_massL = (incoming_lepton_mass > 1e-9) ? incoming_lepton_mass : fl[2].Mass();
            m_massp = Flavour(kf_p_plus).Mass();
            m_massn = Flavour(kf_n).Mass();
            m_massavg = 0.5 * (m_massp+m_massn);

            msg_Out() << "  -> initializing W form factors for incoming nucleon: " << m_nucleon << "\n";
            p_ff_W = std::make_unique<FormFactor_EMnucleon>(incomingboson::W, m_nucleon);
        }

        double operator()(const ATOOLS::Vec4D_Vector &mom) override
        {
            Kinematics k = ComputeKinematics(mom);
            if (k.Q2 < s_Q2min) return 0.0;

            // msg_Out() << " inv mass from mom: m_ki^2" << ki.Abs2() << "m_pi ^ 2 " << pi.Abs2() << "m_kf^2" << kf.Abs2() << "m_pf^2" << pf.Abs2() << std::endl;
            // Mandelstraam variables
            const double s = k.s;
            const double t = k.t; // t < 0 virtual boson exchange
            const double u = k.u; // u very small, so dont use in calculation to avoid rounding errors.

            const double Q2 = k.Q2; // Q2 > 0 for virtual boson exchange
            msg_Tracking()  << METHOD << ": theta = " << k.theta << ", phi = " << k.phi << "\n"
                            << METHOD << ": Q2 = " << Q2 << "\n"
                            << METHOD << ": s = " << s << ", t = " << t << ", u = " << u << "\n"
                            << METHOD << ": ki = " << k.ki << ", pi = " << k.pi << ", kf = " << k.kf << ", pf = " << k.pf << "\n"
                            << METHOD << ": Mass_neutron = " << m_massn << ", Mass_proton = " << m_massp << "\n"
                            << METHOD << ": MassW^2       = " << m_massWsqr << "\n";

            // Get form factors (constants for given Q2)
            NucleonFormFactors ff = p_ff_W->GetFormFactors(Q2);
            const double F1 = ff.F1;
            const double F2 = ff.F2;
            const double FA = ff.FA;
            const double FP = ff.FP;
            msg_Tracking() << METHOD << ": F1 = " << F1 << ", F2 = " << F2 << ", FA = " << FA << ", FP = " << FP << "\n"
                            << METHOD << ": gW^4 = " << m_gW4 << std::endl;

            // precalcs
            const double t2 = t * t;
            const double s2 = s * s;
            const double u2 = u * u;

            // Calc
            const double m  = m_massavg;
            const double me = m_massL;
            const double mp = m_massp;
            const double mn = m_massn;
            
            const double me2 = me * me;
            const double mp2 = mp * mp;
            const double mn2 = mn * mn;
            
            const double calc =
                (me2 - t) * (
                    4.0 * F2 * m * mn * (F1 + FA) * (me2 + mp2 - s)
                  - F2 * (mn2 - s) * (F2 * (mp2 - u) - 4.0 * m * mp * (F1 + FA))
                  - 4.0 * m * mn * (F1 - FA) * (2.0 * m * mp * (F1 + FA) + F2 * (u - mp2))
                  - me2 * (mn * mp * (3.0 * F2 * F2 + FP * FP)
                         + 0.5 * (F2 - FP) * (F2 + FP) * (mn2 + mp2 - t))
                )
              + 2.0 * me2 * (
                    (mn2 - s) * (-2.0 * F2 * m * mp * (F1 + 2.0 * FA)
                               + F2 * F2 * (mp2 - u)
                               - 2.0 * FA * FP * m * mp)
                  - 2.0 * m * mn * (mp2 - u) * (F1 * F2 + FA * (FP - 2.0 * F2))
                )
              + (me2 + mn2 - u) * (
                    4.0 * F2 * m * mp * (F1 - FA) * (me2 - t)
                  + 4.0 * m * m * (F1 - FA) * (F1 - FA) * (mp2 - u)
                  - F2 * F2 * (me2 - t) * (me2 + mp2 - s)
                )
              + 4.0 * m * m * (F1 + FA) * (F1 + FA) * (mn2 - s) * (me2 + mp2 - s)
              + 2.0 * F2 * F2 * mn * mp * (me2 - t) * (me2 - t);
            // Coefficient
            const double prefactor = (4 * m_gW4) / (m_helicityavg * m_massavg*m_massavg * (m_massWsqr + Q2)*(m_massWsqr + Q2));

            double M_squared = prefactor * calc;
            msg_Tracking() << METHOD << ": M2 = " << M_squared << std::endl;

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
    if (!fl[0].IsLepton()  || !fl[1].IsNucleon()) {

        return NULL;
    }
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
        return new EXTRAXS::LN_LN_cc(args, nucleon, false);
    }
    if (fl[0].IsNeutrino() && fl[2].IsChargedLepton()) {
        // Reverse CC: nu N -> l N'  — nucleon identity also reverses
        msg_Out() << "  -> dispatching to LN_LN_cc \n";
        return new EXTRAXS::LN_LN_cc(args, nucleon, true);
    }

    // ---- NC: lepton identity preserved (fl[0]==fl[2], fl[1]==fl[3]) ----
    if (fl[0] != fl[2] || fl[1] != fl[3]) return NULL;

    if (fl[0].IsChargedLepton()) {
        msg_Out() << "  -> dispatching to LN_LN_nc_charged \n";
        return new EXTRAXS::LN_LN_nc_charged(args, nucleon);
    }
    if (fl[0].IsNeutrino()) {
        msg_Out() << "  -> dispatching to LN_LN_nc_neutrino \n";
        return new EXTRAXS::LN_LN_nc_neutrino(args, nucleon);
    }

    return NULL;
}
