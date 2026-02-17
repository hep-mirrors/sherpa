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
    // NOTE: that this class also does nv_pe scattering when the nucleon is a neutron
    //       and that all refernces to protons will become neutrons if that is the incoming nucleon
    // NOTE: this class also works for muon nucleon scattering, as lepton mass is set
    //        by incoming lepton flavour

    class pe_nv : public ME2_Base
    {
    private:
        // Maps physical index to Sherpa's reordered momentum index
        int m_order;
        inline size_t GetMomentumIndex(size_t input_idx) const {
            if (m_order == 0) {
                // proton case: p[0] e-[1] -> ve[2] n[3]
                // ki=1, pi=0, kf=2, pf=3
                static const size_t map[4] = {1,0,2,3};
                return map[input_idx];
            } else if (m_order == 1) {
                // neutron case: ve[0] n[1] -> p[2] e-[3]
                // ki=0, pi=1, kf=3, pf=2
                // Swap initial and final state
                static const size_t map[4] = {0,1,3,2};
                return map[input_idx];
            } else {
                // neutron case: n[0] ve[1] -> p[2] e-[3]
                // ki=0, pi=1, kf=3, pf=2
                // Swap initial and final state
                static const size_t map[4] = {1, 0, 3, 2};
                return map[input_idx];
            }
        }
        
    public:
        pe_nv(const External_ME_Args &args, const incomingnucleon::code &nucleon,
              const Flavour &fl0, const Flavour &fl1, const Flavour &fl2, const Flavour &fl3, int &order);

        double operator()(const ATOOLS::Vec4D_Vector &mom);
        double m_alpha, m_alphas, m_g, m_gA, m_massW;
        Flavour m_fl0, m_fl1, m_fl2, m_fl3;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
        incomingnucleon::code m_nucleon; // neutron or proton decides ordering
    };

    pe_nv::pe_nv(const External_ME_Args &args, const incomingnucleon::code &nucleon,
                 const Flavour &fl0, const Flavour &fl1, const Flavour &fl2, const Flavour &fl3, int &order)
        : ME2_Base(args), m_nucleon(nucleon), m_fl0(fl0), m_fl1(fl1), m_fl2(fl2), m_fl3(fl3), m_order(order)
    {
        msg_Out() << "pe_nv::pe_nv(): Constructor called" << std::endl;
        msg_Out() << "  with flavours" << std::endl;
        msg_Out() << "  fl[0]:" << m_fl0 << " fl[1]:" << m_fl1 << " fl[2]:" << m_fl2 << " fl[3]:" << m_fl3 << std::endl;
        msg_Out() << "  m_order = " << m_order << std::endl;

        m_sintt = 2;  // T-channel only , needed for nv-pe to integrate...
        // m_sintt = 0; // leads to no channels?
        m_oew = 0;   // EW order zero (no loops)
        m_oqcd = 0;  // QCD order zero
        msg_Out() << "pe_nv::pe_nv(): Getting constants" << std::endl;
        m_alpha = 1./137; // m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");
        // m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");
        m_massW = Flavour(kf_Wplus).Mass();
        m_g = 0.653;            // MODEL::s_model->ScalarConstant("weak_cpl");
        m_gA = 1.27; //axial coupling

        msg_Out() << "pe_nv::pe_nv(): Creating FormFactor_EMnucleon for W" << std::endl;
        p_formfactor = std::make_unique<FormFactor_EMnucleon>(incomingboson::W, m_nucleon);

        msg_Out() << "pe_nv::pe_nv(): Constructor complete for nucleon " << m_nucleon << std::endl;
    }

    double pe_nv::operator()(const ATOOLS::Vec4D_Vector &momenta)
    {
        // SHERPA ORDER indices
        // msg_Out()<<"pe_nv::operator(): Called with " << std::endl;
        // k: lepton, p: hadron
        // Use GetMomentumIndex to map input indices to Sherpa's ordering
        const auto &ki = momenta[GetMomentumIndex(0)];  // initial lepton
        const auto &pi = momenta[GetMomentumIndex(1)];  // initial hadron
        const auto &kf = momenta[GetMomentumIndex(2)];  // final lepton
        const auto &pf = momenta[GetMomentumIndex(3)];  // final hadron

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

        msg_Out() << "pe_nv::operator(): theta = " << theta << std::endl;
        msg_Out() << "pe_nv::operator(): phi = " << phi << std::endl;

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
        const double s = (ki + pi).Abs2();
        const double t = (ki - kf).Abs2(); // t < 0 virtual boson exchange
        const double u = (ki - pf).Abs2(); // u very small, so dont use in calculation to avoid rounding errors.

        const double Q2 = -t; // Q2 > 0 for virtual boson exchange
        msg_Out() << "pe_nv::operator(): Q2 = " << Q2 << std::endl;
        msg_Out() << "pe_nv::operator(): s = " << s << ", t = " << t << ", u = " << u << std::endl;

        // Add a small cutoff to avoid IR divergence at Q2=0 (forward scattering)
        // cutoff order approx mass of electron, we dont expect slow moving electrons
        if (Q2 < 1e-5)
        {
            // msg_Out() << "pe_nv::operator(): Q2 < Q2_min, returning 0.0" << std::endl;
            return 0.0;
        }

        // Get form factors (constants for given Q2)
        NucleonFormFactors ff = p_formfactor->GetFormFactors(Q2);
        const double F1 = ff.F1;
        const double F2 = ff.F2;
        const double FA = ff.FA;
        const double FP = ff.FP;
        msg_Out() << "pe_nv::operator(): F1 = " << F1 << ", F2 = " << F2 << ", FA = "<< FA << ", FP" << FP << std::endl;

        // precalcs
        const double MW2 = m_massW * m_massW;
        const double gW4 = (m_g*m_g*m_g*m_g)/4;
        const double sum_mass = Mn+Mp;
        const double sum_mass2 = (Mn + Mp)*(Mn+Mp);
        const double diff_mass = Mn-Mp;
        const double diff_mass2 = (Mn-Mp)*(Mn-Mp);
        const double t2 = t*t;
        const double s2 = s*s;
        const double st = s*t;

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

        msg_Out() << "pe_nv::operator(): M_squared = " << M_squared << "\n" << std::endl;
        return M_squared;
    }
}

DECLARE_TREEME2_GETTER(EXTRAXS::pe_nv, "pe_nv")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::pe_nv>::
operator()(const External_ME_Args &args) const
{
    const Flavour_Vector &fl = args.Flavours();

    // check if 2-2 CC scattering
    if (fl.size() != 4)
        return NULL;
    // check lepton family (any ordering)
    if (fl[1].LeptonFamily() != fl[2].LeptonFamily()) {
        if (fl[1].LeptonFamily() != fl[3].LeptonFamily())
        {
            return NULL;
        }
    }

    // p e- -> n ve (after Sherpa reordering: p e- -> ve n)
    if (fl[0] == Flavour(kf_p_plus) && fl[1].IsChargedLepton() && 
         (fl[3] == Flavour(kf_n) && fl[2].IsNeutrino()))
    {
        int m_order = 0;
        return new pe_nv(args, incomingnucleon::proton, fl[0], fl[1], fl[2], fl[3], m_order);
    } 
    // n ve -> p e- (after Sherpa reordering: ve n  -> p e- )
    else if (fl[1] == Flavour(kf_n) && fl[0].IsNeutrino() && 
             (fl[2] == Flavour(kf_p_plus) && fl[3].IsChargedLepton()))
    {
        int m_order = 1;
        return new pe_nv(args, incomingnucleon::neutron, fl[0], fl[1], fl[2], fl[3], m_order);
    }
    else if (fl[0] == Flavour(kf_n) && fl[1].IsNeutrino() &&
             (fl[2] == Flavour(kf_p_plus) && fl[3].IsChargedLepton()))
    {
        int m_order = 2;
        return new pe_nv(args, incomingnucleon::neutron, fl[0], fl[1], fl[2], fl[3], m_order);
    }
    else
    {
        return NULL; // Not an appropriate process
    }

    return NULL;
}
