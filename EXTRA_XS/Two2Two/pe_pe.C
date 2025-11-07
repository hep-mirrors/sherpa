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

    class pe_pe : public ME2_Base
    {
    private:

    public:
        pe_pe(const External_ME_Args &args, const incomingboson::code &boson, const incomingnucleon::code &nucleon);

        double operator()(const ATOOLS::Vec4D_Vector &mom);
        double m_alpha, m_alphas, m_s, Me2, Mp2;
        Flavour m_flv;
        std::unique_ptr<FormFactor_EMnucleon> p_formfactor;
        incomingboson::code m_boson;
        incomingnucleon::code m_nucleon;
    };

    pe_pe::pe_pe(const External_ME_Args &args, const incomingboson::code &boson, const incomingnucleon::code &nucleon)
        : ME2_Base(args), m_boson(boson), m_nucleon(nucleon)
    {
        msg_Out()<<"pe_pe::pe_pe(): Constructor called"<<std::endl;
        m_sintt = 2;  // T-channel only , maybe look at single channel restriction later
        //m_sintt = 0;  // No s/t/u channel restriction
        m_oew = 0; // EW order zero (no loops)
        m_oqcd = 0; // QCD order zero 
        msg_Out()<<"pe_pe::pe_pe(): Getting alpha_QED"<<std::endl;
        m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");
        msg_Out()<<"pe_pe::pe_pe(): Getting strong_cpl"<<std::endl;
        m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");

        msg_Out()<<"pe_pe::pe_pe(): Creating FormFactor_EMnucleon"<<std::endl;
        p_formfactor = std::make_unique<FormFactor_EMnucleon>(m_boson, m_nucleon);
        msg_Out()<<"pe_pe::pe_pe(): Constructor complete"<<std::endl;
        // m_flvs = args.Flavours();
        // if (m_flvs[0].Charge() == m_flvs[1].Charge())
        //     m_ss = true;
        // else
        //     m_ss = false;
        msg_Out()<<"pe_pe::pe_pe(): Setting Me2 and Mp2"<<std::endl;
        Me2 = pow(Flavour(kf_e).Mass(), 2);
        Mp2 = pow(Flavour(kf_p_plus).Mass(), 2);
    }

    double pe_pe::operator()(const ATOOLS::Vec4D_Vector &momenta)
    {
        // indices: # p[0] e-[1] -> p[2] e-[3]
        // k: lepton, p: hadron
        const auto &ki = momenta[1];
        const auto &pi = momenta[0];
        const auto &kf = momenta[3];
        const auto &pf = momenta[2];
        
        // Mandelstraam variables
        double s = (ki + pi).Abs2();
        double t = (ki - kf).Abs2();  // t < 0 for spacelike
        double u = (ki - pf).Abs2();
        
        double Q2 = -t;  // Q2 > 0 for spacelike photon exchange
        
        // Add a small cutoff to avoid IR divergence at Q2=0 (forward scattering)
        // For low-energy fixed target, use a smaller cutoff
        const double Q2_min = 1e-10;
        if (Q2 < Q2_min) {
            msg_Out() << "pe_pe::operator(): Q2 < Q2_min, returning 0.0" << std::endl;
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
        double pre_factor = 32 * M_PI * M_PI * sqr((*aqed)(m_s)) / (Q2 * Q2);

        return pre_factor * L_munu_H_munu;
    }
}

DECLARE_TREEME2_GETTER(EXTRAXS::pe_pe, "pe_pe")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base, PHASIC::External_ME_Args, EXTRAXS::pe_pe>::
operator()(const External_ME_Args &args) const
{
    const Flavour_Vector fl = args.Flavours();
    if (fl.size() != 4) return NULL;
    
    // Check for P+ e- -> P+ e- : sherpa orders hadrons first!!!!!
    // Sherpas initial state: fl[0]=P+, fl[1]=e-
    // Sherpas final state: fl[2]=P+, fl[3]=e-
    if (fl[0] == Flavour(kf_p_plus) && fl[1] == Flavour(kf_e) &&
        fl[2] == Flavour(kf_p_plus) && fl[3] == Flavour(kf_e))
    {
        return new pe_pe(args, incomingboson::photon, incomingnucleon::proton);
    }
    return NULL;
}
