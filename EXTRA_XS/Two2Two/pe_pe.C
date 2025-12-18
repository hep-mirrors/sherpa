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
              const Flavour &lepton_fl, const Flavour &nucleon_fl, 
              bool include_Z_interference = true);
        
        double operator()(const ATOOLS::Vec4D_Vector &mom);
        double m_alpha, m_alphas, m_e2, m_massZ, m_massW; 
        Flavour m_fl_lepton, m_fl_nucleon;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor_Z;  
        incomingnucleon::code m_nucleon;
        bool m_include_Z;  // Include photon-Z interference for NC
    };

    pe_pe::pe_pe(const External_ME_Args &args, const incomingnucleon::code &nucleon,
                 const Flavour &lepton_fl, const Flavour &nucleon_fl,
                 bool include_Z_interference)
        : ME2_Base(args), m_nucleon(nucleon), m_fl_lepton(lepton_fl), 
          m_fl_nucleon(nucleon_fl), m_include_Z(include_Z_interference)
    {
        msg_Out()<<"pe_pe::pe_pe(): Constructor called"<<std::endl;
        m_sintt = 2;  // T-channel only , maybe look at single channel restriction later
        //m_sintt = 0;  // No s/t/u channel restriction
        m_oew = 0; // EW order zero (no loops)
        m_oqcd = 0; // QCD order zero 
        msg_Out()<<"pe_pe::pe_pe(): Getting constants"<<std::endl;
        m_alpha = 1.0/137.0; //MODEL::s_model->ScalarConstant("fine_structure_const");
        m_e2 = 4 * M_PI * m_alpha;
        //m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");
        m_massZ = Flavour(kf_Z).Mass();

        msg_Out()<<"pe_pe::pe_pe(): Creating FormFactor_EMnucleon for photon"<<std::endl;
        p_formfactor = std::make_unique<FormFactor_EMnucleon>(incomingboson::photon, m_nucleon);
        
        // If including Z interference, create Z form factor too
        if (m_include_Z) {
            msg_Out()<<"pe_pe::pe_pe(): Also creating FormFactor_EMnucleon for Z (interference)"<<std::endl;
            p_formfactor_Z = std::make_unique<FormFactor_EMnucleon>(incomingboson::Z, m_nucleon);
        }
        
        msg_Out()<<"pe_pe::pe_pe(): Constructor complete for nucleon "<< m_nucleon <<std::endl;
    }

    double pe_pe::operator()(const ATOOLS::Vec4D_Vector &momenta)
    {
        msg_Out()<<"pe_pe::operator(): Called with " << std::endl;
        // indices: # p[0] e-[1] -> p[2] e-[3]
        // k: lepton, p: hadron
        const auto &ki = momenta[1];
        const auto &pi = momenta[0];
        const auto &kf = momenta[3];
        const auto &pf = momenta[2];

        msg_Out()<<"pe_pe::operator(): ki = " << ki << std::endl;
        msg_Out()<<"pe_pe::operator(): pi = " << pi << std::endl;
        msg_Out()<<"pe_pe::operator(): kf = " << kf << std::endl;
        msg_Out()<<"pe_pe::operator(): pf = " << pf << std::endl;

        //masses squared (on shell)
        const double Me2 = m_fl_lepton.Mass() * m_fl_lepton.Mass();
        const double Mp2 = m_fl_nucleon.Mass() * m_fl_nucleon.Mass();

        msg_Out()<<"pe_pe::operator(): Me2 = " << Me2 << std::endl;
        msg_Out()<<"pe_pe::operator(): Mp2 = " << Mp2 << std::endl; 
        
        // Mandelstraam variables
        double s = (ki + pi).Abs2();
        double t = (ki - kf).Abs2();  // t < 0 virtual boson exchange
        double u = (ki - pf).Abs2();

        msg_Out()<<"pe_pe::operator(): s = " << s << std::endl;
        msg_Out()<<"pe_pe::operator(): t = " << t << std::endl;
        msg_Out()<<"pe_pe::operator(): u = " << u << std::endl;
        
        double Q2 = -t;  // Q2 > 0 for virtual boson exchange

        // Add a small cutoff to avoid IR divergence at Q2=0 (forward scattering)
        // cutoff order approx mass of electron, we dont expect slow moving electrons
        const double Q2_min = 1e-3;
        if (Q2 < Q2_min) {
            // msg_Out() << "pe_pe::operator(): Q2 < Q2_min, returning 0.0" << std::endl;
            return 0.0;
        }

        NucleonFormFactors ff = p_formfactor->GetFormFactors(Q2);
        double F1 = ff.F1;
        double F2 = ff.F2;
        double F12 = F1 + F2;

        double A = (s - Me2 - Mp2)/2.0;  // (ki . pi)
        double B = -(u - Me2 - Mp2)/2.0; // (ki . pf)
        double C = Me2 - (t)/2.0;       // (ki . kf)
        double D = Mp2 - (t)/2.0;       // (pi . pf)

        // Calculate spin averaged matrix element squared
        double term1coeff = 4 * F12 * F12;
        double term1 = term1coeff * (2 * A * A + 2 * B * B + C * t + D * t);
        double term2coeff = 2*( (F2*F2* (Mp2 + 4*D))/(8*Mp2) -2*F2*F12 );
        double term2 = term2coeff * (2*( A*A + B*B +2*A*B) + (C + D +2*B)*t*0.5);

        double L_munu_H_munu = - (term1 + term2);
        double pre_factor = 2* m_e2 * m_e2 / (Q2 * Q2);

        double M_squared = pre_factor * L_munu_H_munu;

        if (m_include_Z && p_formfactor_Z)
        {
            NucleonFormFactors ff_Z = p_formfactor_Z->GetFormFactors(Q2);
            double F1_Z = ff_Z.F1;
            double F2_Z = ff_Z.F2;
            double FA_Z = ff_Z.FA;
            double FP_Z = ff_Z.FP;
            double F12_Z = F1_Z + F2_Z;

            double A_Z = (s - Me2 - Mp2) / 2.0;  // (ki . pi)
            double B_Z = -(u - Me2 - Mp2) / 2.0; // (ki . pf)
            double C_Z = Me2 - (t) / 2.0;        // (ki . kf)
            double D_Z = Mp2 - (t) / 2.0;        // (pi . pf)

            // Calculate spin averaged matrix element squared
            double term1coeff_Z = 4 * F12_Z * F12_Z;
            double term1_Z = term1coeff_Z * (2 * A_Z * A_Z + 2 * B_Z * B_Z + C_Z * t + D_Z * t);
            double term2coeff_Z = 2 * ((F2_Z * F2_Z * (Mp2 + 4 * D_Z)) / (8 * Mp2) - 2 * F2_Z * F12_Z);
            double term2_Z = term2coeff_Z * (2 * (A_Z * A_Z + B_Z * B_Z + 2 * A_Z * B_Z) + (C_Z + D_Z + 2 * B_Z) * t * 0.5);
            // plus axial terms 

            double L_munu_H_munu_Z = -(term1_Z + term2_Z);
            double pre_factor_Z = 2 * m_e2 * m_e2 / pow((Q2+m_massZ*m_massZ),2);

            double M_squared_Z = pre_factor_Z * L_munu_H_munu_Z;

            // interference term 
            double interference_term = 0.0; // to be calculated

            // all together 
            M_squared += M_squared_Z + interference_term;
        }

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

