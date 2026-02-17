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

namespace EXTRAXS {
    // NOTE: that this class also does ne-ne scattering when the nucleon is a neutron
    //       and that all refernces to protons will become neutrons if that is the incoming nucleon
    // NOTE: this class also works for muon nucleon scattering, as lepton mass is set 
    //        by incoming lepton flavour

    class pe_pe : public ME2_Base
    {
    private:

    public:
        pe_pe(const External_ME_Args &args, const incomingnucleon::code &nucleon, 
              const Flavour &lepton_fl, const Flavour &nucleon_fl);
        
        double operator()(const ATOOLS::Vec4D_Vector &mom);
        double m_alpha, m_alphas, m_e2, m_massZ; 
        Flavour m_fl_lepton, m_fl_nucleon;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor_Z;  
        incomingnucleon::code m_nucleon;
        bool m_include_Z;  // Include photon-Z interference for NC
    };

    pe_pe::pe_pe(const External_ME_Args &args, const incomingnucleon::code &nucleon,
                 const Flavour &lepton_fl, const Flavour &nucleon_fl)
        : ME2_Base(args), m_nucleon(nucleon), m_fl_lepton(lepton_fl), 
          m_fl_nucleon(nucleon_fl)
    {
        msg_Out()<<"pe_pe::pe_pe(): Constructor called"<<std::endl;
        m_sintt = 2;  // T-channel only , maybe look at single channel restriction later
        //m_sintt = 0;  // No s/t/u channel restriction
        m_oew = 0; // EW order zero (no loops)
        m_oqcd = 0; // QCD order zero 
        // msg_Out()<<"pe_pe::pe_pe(): Getting constants"<<std::endl;
        m_alpha = 1.0/137.0; //MODEL::s_model->ScalarConstant("fine_structure_const");
        m_e2 = 4 * M_PI * m_alpha;
        //m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");
        m_massZ = Flavour(kf_Z).Mass();
        
        // Read Include_Z_Interference from yaml settings 
        ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
        m_include_Z = s["Form_Factor"]["Include_Z_Interference"].SetDefault(true).Get<bool>();

        //TODO: put settings in getter instead of hardcoding

        msg_Out()<<"pe_pe::pe_pe(): Creating FormFactor_EMnucleon for photon"<<std::endl;
        p_formfactor = std::make_unique<FormFactor_EMnucleon>(incomingboson::photon, m_nucleon);
        
        // If including Z interference, create Z form factor too
        if (m_include_Z) {
            msg_Out()<<"pe_pe::pe_pe(): Also creating FormFactor_EMnucleon for Z (interference)"<<std::endl;
            p_formfactor_Z = std::make_unique<FormFactor_EMnucleon>(incomingboson::Z, m_nucleon);
        }
        
        msg_Out()<<"pe_pe::pe_pe(): Constructor complete for nucleon: "<<m_nucleon<<std::endl;
    }

    double pe_pe::operator()(const ATOOLS::Vec4D_Vector &momenta)
    {
        // msg_Out()<<"pe_pe::operator(): Called with " << std::endl;
        // indices: # p[0] e-[1] -> p[2] e-[3]
        // k: lepton, p: hadron
        const auto &ki = momenta[1];
        const auto &pi = momenta[0];
        const auto &kf = momenta[3];
        const auto &pf = momenta[2];

        // msg_Out() << "pe_pe::operator(): ki = " << ki << std::endl;
        // msg_Out() << "pe_pe::operator(): pi = " << pi << std::endl;
        // msg_Out() << "pe_pe::operator(): kf = " << kf << std::endl;
        // msg_Out() << "pe_pe::operator(): pf = " << pf << std::endl;
        
        // Calculate scattering angle theta between initial and final lepton
        Vec3D ki_vec(ki);
        Vec3D kf_vec(kf);
        double cos_theta = (ki_vec * kf_vec) / (ki_vec.Abs() * kf_vec.Abs());
        double theta = acos(cos_theta);
        double phi = atan2(kf[2], kf[1]);
        
        // msg_Out() << "pe_pe::operator(): theta = " << theta << std::endl;
        // msg_Out() << "pe_pe::operator(): phi = " << phi << std::endl;

        //masses squared (on shell)
        const double Me2 = m_fl_lepton.Mass() * m_fl_lepton.Mass();
        const double Mp2 = m_fl_nucleon.Mass() * m_fl_nucleon.Mass();
        // msg_Out() << "pe_pe::operator(): Me2 = " << Me2 << ", Mp2 = " << Mp2 << std::endl;
        
        // Mandelstraam variables
        const double s = (ki + pi).Abs2();
        const double t = (ki - kf).Abs2();  // t < 0 virtual boson exchange
        const double u = (ki - pf).Abs2();
        
        const double Q2 = -t;  // Q2 > 0 for virtual boson exchange
        // msg_Out() << "pe_pe::operator(): Q2 = " << Q2 << std::endl;
        // msg_Out() << "pe_pe::operator(): s = " << s << ", t = " << t << ", u = " << u << std::endl;

        // Add a small cutoff to avoid IR divergence at Q2=0 (forward scattering)
        // cutoff order approx mass of electron, we dont expect slow moving electrons
        if (Q2 < 1e-5)
        {
            // msg_Out() << "pe_pe::operator(): Q2 < Q2_min, returning 0.0" << std::endl;
            return 0.0;
        }

        // Get form factors (constants for given Q2)
        NucleonFormFactors ff = p_formfactor->GetFormFactors(Q2);
        const double F1 = ff.F1;
        const double F2 = ff.F2;
        // msg_Out() << "pe_pe::operator(): F1 = " << F1 << ", F2 = " << F2 << std::endl;
        
        // Pre-compute common combinations
        const double s2 = s*s;
        const double u2 = u*u;
        const double t2 = t*t;
        const double Me4 = Me2 * Me2;
        const double Me6 = Me4 * Me2;
        const double Mp4 = Mp2 * Mp2;
        const double F1_2 = F1 * F1;
        const double F2_2 = F2 * F2;
        const double twoMe2 = 2*Me2;
        const double twoMp2 = 2*Mp2;
        const double twoMe2_plus_t = twoMe2 + t;
        const double twoMe2_plus_twoMp2_minus_s_minus_u = twoMe2 + twoMp2 - s - u;

        double calc = 
            2*F1_2*Mp2 * (
                2*Me4
                + twoMe2*(twoMp2 - s + t - u)
                + 2*Mp4
                - twoMp2*(s - t + u)
                + s2 + u2
            )
            + 4*F1*F2*Mp2*twoMe2_plus_t*twoMe2_plus_twoMp2_minus_s_minus_u
            + F2_2 * (
                4*Me6
                + Me4*(8*Mp2 - 4*s + t - 4*u)
                + Me2*(4*Mp4
                       + Mp2*(-4*s + 6*t - 4*u)
                       + s2 - s*t + 2*s*u - t2 - t*u + u2)
                + t*(Mp4 - Mp2*(s - t + u) + s*u)
            );

        const double Q4 = Q2 * Q2;
        const double pre_factor = m_e2 * m_e2 / (Mp2 * Q4);
        double M_squared = pre_factor * calc;
        // msg_Out() << "pe_pe::operator(): M_squared (photon only) = " << M_squared << std::endl;

        if (m_include_Z)
        {
            NucleonFormFactors ff_Z = p_formfactor_Z->GetFormFactors(Q2);
            const double F1_Z = ff_Z.F1;
            const double F2_Z = ff_Z.F2;
            const double FA_Z = ff_Z.FA;
            const double FP_Z = ff_Z.FP;
            // msg_Out() << "pe_pe::operator(): F1_Z = " << F1_Z << ", F2_Z = " << F2_Z 
            //           << ", FA_Z = " << FA_Z << ", FP_Z = " << FP_Z << std::endl;

            // Pre-compute common Z combinations
            const double Mp6 = Mp4 * Mp2;
            const double sW2 = p_formfactor_Z->m_sin2thetaW;
            const double cW2 = p_formfactor_Z->m_cos2thetaW;
            const double cW4 = cW2 * cW2;
            const double g2 = p_formfactor_Z->m_g * p_formfactor_Z->m_g;
            const double g4 = g2 * g2;
            const double gVe = -0.5 + 2*sW2;
            const double gAe = -0.5;
            const double gVe2 = gVe * gVe;
            const double gAe2 = gAe * gAe;
            const double FAZ2 = FA_Z * FA_Z;
            const double F1Z2 = F1_Z * F1_Z;
            const double F2Z2 = F2_Z * F2_Z;
            const double FPZ2 = FP_Z * FP_Z;
            const double MZ2 = m_massZ * m_massZ;
            const double MZ2_plus_Q2 = MZ2 + Q2;
            const double MZ2_plus_Q2_sq = MZ2_plus_Q2 * MZ2_plus_Q2;
            
            // Compute interference term numerator (simplified with constants)
            const double s_minus_u = s - u;
            const double gAe2_plus_gVe2 = gAe2 + gVe2;
            
            double calc_Z = 
                - 2 * cW2 * g2 * Q2 * MZ2_plus_Q2 * (-1.0) * (  // Qf = -1.0 for electron
                    2 * F1 * (
                        twoMe2_plus_twoMp2_minus_s_minus_u *
                            (F2_Z * gVe * twoMe2_plus_t + FA_Z * gAe * (u - s))
                        + F1_Z * gVe * (
                            2*Me4
                            + twoMe2*(twoMp2 - s + t - u)
                            + 2*Mp4
                            + s2 + u2
                            - twoMp2*(s - t + u)
                        )
                    ) * Mp2

                    + F2 * (
                        2 * (
                            F1_Z*gVe*twoMe2_plus_t*twoMe2_plus_twoMp2_minus_s_minus_u
                            + FA_Z*gAe*t*s_minus_u
                        ) * Mp2

                        + F2_Z * gVe * (
                            4*Me6
                            + (8*Mp2 - 4*s + t - 4*u)*Me4
                            + (4*Mp4 + (-4*s + 6*t - 4*u)*Mp2 + s2 - t2 + u2 - s*t + 2*s*u - t*u)*Me2
                            + t*(Mp4 - (s - t + u)*Mp2 + s*u)
                        )
                    )
                ) * m_e2

                + g4 * Q4 * (
                    4*FAZ2*gAe2_plus_gVe2*Mp6
                    + 40*FAZ2*gAe2*twoMe2*Mp4
                    + 16*FA_Z*FP_Z*gAe2*twoMe2*Mp4
                    - 8*FAZ2*gVe2*twoMe2*Mp4
                    - 4*FAZ2*gAe2_plus_gVe2*s*Mp4
                    - 4*FAZ2*gAe2_plus_gVe2*t*Mp4

                    + 4*FAZ2*gAe2*Me4*Mp2
                    + 16*FA_Z*FP_Z*gAe2*Me4*Mp2
                    + 4*FAZ2*gVe2*Me4*Mp2

                    + 2*FAZ2*gAe2_plus_gVe2*s2*Mp2
                    + 2*FAZ2*gAe2_plus_gVe2*u2*Mp2

                    - 4*FAZ2*gAe2*twoMe2*s*Mp2
                    - 8*FA_Z*FP_Z*gAe2*twoMe2*s*Mp2
                    - 4*FAZ2*gVe2*twoMe2*s*Mp2

                    - 4*FAZ2*gAe2*twoMe2*t*Mp2
                    + 4*FAZ2*gVe2*twoMe2*t*Mp2

                    - 4*F2_Z * (
                        F1_Z * (gAe2*(4*Me2 - t) - gVe2*twoMe2_plus_t) * twoMe2_plus_twoMp2_minus_s_minus_u
                        + 2*FA_Z*gAe*gVe * t * s_minus_u
                    ) * Mp2

                    - 8*F1_Z*FA_Z*gAe*gVe * twoMe2_plus_twoMp2_minus_s_minus_u * s_minus_u * Mp2

                    - 4*FA_Z*(2*FP_Z*gAe2*Me2 + FA_Z*gAe2_plus_gVe2*(Me2 + Mp2)) * u * Mp2

                    + 2*F1Z2 * (
                        (2*Me4 - twoMe2*(twoMp2 + s + t + u) + 2*Mp4 + s2 + u2 - twoMp2*(s - t + u)) * gAe2
                        + gVe2 * (
                            2*Me4 + twoMe2*(twoMp2 - s + t - u)
                            + 2*Mp4 + s2 + u2
                            - twoMp2*(s - t + u)
                        )
                    ) * Mp2

                    + FPZ2 * gAe2 * Me2 * t2

                    - F2Z2 * (
                        gAe2 * (
                            4*Me6
                            + (8*Mp2 - 4*s - t - 4*u)*Me4
                            + (4*Mp4 + 2*(-2*s + t - 2*u)*Mp2 + (s + u)*(s + t + u))*Me2
                            + t * (-Mp4 + (s - t + u)*Mp2 - s*u)
                        )
                        - gVe2 * (
                            4*Me6
                            + (8*Mp2 - 4*s + t - 4*u)*Me4
                            + (4*Mp4 + (-4*s + 6*t - 4*u)*Mp2 + s2 - t2 + u2 - s*t + 2*s*u - t*u)*Me2
                            + t*(Mp4 - (s - t + u)*Mp2 + s*u)
                        )
                    )
                );

            // full 2*M*M_Z + M2_Z interference terms 
            double interference_term = calc_Z / (cW4 * Mp2 * Q4 * MZ2_plus_Q2_sq); 

            // all together 
            M_squared +=  interference_term;
        //     msg_Out() << "pe_pe::operator(): M_squared (with Z interference) = " << M_squared << std::endl;
        }
        // msg_Out() << "pe_pe::operator(): M_squared = " << M_squared << std::endl;
        return M_squared;
    }
}

DECLARE_TREEME2_GETTER(EXTRAXS::pe_pe, "pe_pe")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::pe_pe>::
operator()(const External_ME_Args &args) const
{
    const Flavour_Vector &fl = args.Flavours();
    incomingnucleon::code nucleon = incomingnucleon::off;

    // check if elastic NC scattering
    if (fl.size() != 4) return NULL;
    if (fl[0] != fl[2])
        return NULL;
    if (fl[1] != fl[3])
        return NULL;

    // Sherpa orders: fl[0]=hadron_in, fl[1]=lepton_in, fl[2]=hadron_out, fl[3]=lepton_out
    if (fl[0] == Flavour(kf_p_plus)) {
        nucleon = incomingnucleon::proton;
    } else if (fl[0] == Flavour(kf_n)) {
        nucleon = incomingnucleon::neutron;
    } else {
        return NULL; // Not a nucleon
    }
    
    // check leptons for NC interaction (electron, muon, or tau)
    if (fl[1].IsChargedLepton()) {
        return new pe_pe(args, nucleon, fl[1], fl[0]);
    }

    return NULL;
}

