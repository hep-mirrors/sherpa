#include "HADRON_RESCATTERING/Main/Deuteron_Production.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <cmath>

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Deuteron_Coalescence::Deuteron_Coalescence(double b_coal,
                                            double p_coal,
                                            CoalescenceModel model)
  :m_b_coal(b_coal),
    m_p_coal(p_coal),
    m_model(model)
  {
  }


Deuteron_Coalescence::FormationResult Deuteron_Coalescence::CheckPair(Particle* A, Particle* B) {
    FormationResult res;
    Flavour flA = A->Flav(), flB = B->Flav();
    long unsigned int kfA = flA.Kfcode(), kfB = flB.Kfcode();
    bool antiA = flA.IsAnti(), antiB = flB.IsAnti();

    bool is_pn = (kfA==2212 && kfB==2112) || (kfA==2112 && kfB==2212);
    bool is_pp = (kfA==2212 && kfB==2212);
    bool is_nn = (kfA==2112 && kfB==2112);

    if ((!is_pn && !is_pp && !is_nn) || (antiA != antiB)) {
        msg_Out() << "Skipping: not a valid nucleon pair\n";
        return res;
    }

    res.isAnti = antiA;
    res.part1  = A;
    res.part2  = B;

    Vec4D P    = A->Momentum() + B->Momentum();
    res.sqrtS  = P.Abs();
    res.P_vec  = P;

    Poincare boost(P);
    Vec4D p1cm = A->Momentum();
    boost.Boost(p1cm);
    res.p_rel  = sqrt(p1cm[1]*p1cm[1]
                    + p1cm[2]*p1cm[2]
                    + p1cm[3]*p1cm[3]);
    return res;
}

Deuteron_Coalescence::FormationResult Deuteron_Coalescence::EvaluateFormation(Deuteron_Coalescence::FormationResult res,const double& b2_fm2) {
    double Prob = 0.0;

    if (m_model == CoalescenceModel::HardDisk) {
        Prob = (b2_fm2 < m_b_coal*m_b_coal
             && res.p_rel < m_p_coal) ? 1.0 : 0.0;
        res.process = "radiative"; // Previous models always assumed after collision. there was a deuteron and a photon.
    }
    else if (m_model == CoalescenceModel::Exponential) {
        double C_b  = (m_b_coal*m_b_coal) / (0.75*M_PI);
        double C_p  = (m_p_coal*m_p_coal) / (0.75*M_PI);
        Prob        = 0.75 * std::exp(-b2_fm2/C_b)
                          * std::exp(-(res.p_rel*res.p_rel)/C_p);
        res.process = "radiative";//same as harddisk.
    }
    else { // CoalescenceModel::CrossSectionModel case. 
        // get all cross sections
        double sig_rad = m_bb.pnToDGamma(res.part1, res.part2);
        double sig_spi = 0.0, sig_dpi = 0.0;//read it as single pion double pion

        long unsigned int kfA = res.part1->Flav().Kfcode();
        long unsigned int kfB = res.part2->Flav().Kfcode();
        bool is_pn = (kfA==2212 && kfB==2112)||(kfA==2112 && kfB==2212);
        bool is_pp = (kfA==2212 && kfB==2212);
        bool is_nn = (kfA==2112 && kfB==2112);

        if (is_pn) {
            sig_rad = m_bb.pnToDGamma(res.part1, res.part2);
            sig_spi = m_bb.pnToDPiZero(res.part1, res.part2);
            sig_dpi = m_bb.pnToPiZeroPiZero(res.part1, res.part2)
                    + m_bb.pnToPiPlusPiMinus(res.part1, res.part2);
        }
        else if (is_pp) {
            sig_rad = 0.0;
            sig_spi = m_bb.ppToDPiPlus(res.part1, res.part2);
            sig_dpi = m_bb.ppToPiPlusPiZero(res.part1, res.part2);
        }
        else if (is_nn) {
            sig_rad = 0.0;
            sig_spi = m_bb.nnToDPiMinus(res.part1, res.part2);
            sig_dpi = m_bb.nnToPiMinusPiZero(res.part1, res.part2);
        }

        double sig_total = sig_rad + sig_spi + sig_dpi;//this is in microbarns.
        double sigma0    = (1/1.8)*1e6; //taken from table 8 of the paper, which i am not sure how to make this work for different E scale?
        Prob             = std::min(sig_total / sigma0, 1.0);

        // select a process weighted by cross section
        double r = ran->Get() * sig_total;
        if      (r < sig_rad)            res.process = "radiative";
        else if (r < sig_rad + sig_spi)  res.process = "single_pion";
        else                             res.process = "double_pion";

        msg_Out() << "[CrossSection]"
                  << " sig_rad="  << sig_rad
                  << " sig_spi="  << sig_spi
                  << " sig_dpi="  << sig_dpi
                  << " P="        << Prob
                  << " process="  << res.process << "\n";
    }

    msg_Out() << "Candidate: "  << (res.isAnti ? "anti" : "")
              << " p_rel="  << res.p_rel
              << " sqrtS="  << res.sqrtS
              << " Prob="   << Prob << "\n";

    if (Prob <= ran->Get()) return res; // not formed
    res.formed = true;
    return res;
}

Blob* Deuteron_Coalescence::BuildBlob(const FormationResult& res) {
    double m_d = 1.875613;
    if (res.sqrtS < m_d) {
        msg_Out() << "sqrtS=" << res.sqrtS << " < m_d, skipping\n";
        return nullptr;
    }

    Blob* blob = new Blob();
    blob->SetPosition(res.part1->Position());
    res.part1->SetStatus(part_status::decayed);
    res.part2->SetStatus(part_status::decayed);
    blob->AddToInParticles(res.part1);
    blob->AddToInParticles(res.part2);

    double s   = res.sqrtS * res.sqrtS;
    Poincare boost(res.P_vec);

    // isotropic direction
    double cosT = 1.0 - 2.0*ran->Get();
    double sinT = sqrt(std::max(0.0, 1.0-cosT*cosT));
    double phi  = 2.0*M_PI*ran->Get();

    Vec4D p_d_cm, p_X_cm;
    Flavour fl_X;

    if (res.process == "radiative") {
        // 2-body: d + gamma
        double E_d = (s + m_d*m_d) / (2.0*res.sqrtS);
        double E_g = (s - m_d*m_d) / (2.0*res.sqrtS);
        double p_d = sqrt(std::max(E_d*E_d - m_d*m_d, 0.0));
        p_d_cm = Vec4D(E_d,  p_d*sinT*cos(phi),  p_d*sinT*sin(phi),  p_d*cosT);
        p_X_cm = Vec4D(E_g, -p_d*sinT*cos(phi), -p_d*sinT*sin(phi), -p_d*cosT);
        fl_X   = Flavour(kf_photon);
    }
    else if (res.process == "single_pion") {
        // 2-body: d + pi
        double m_pi = 0.139570;
        double E_d  = (s + m_d*m_d - m_pi*m_pi) / (2.0*res.sqrtS);
        double E_pi = (s - m_d*m_d + m_pi*m_pi) / (2.0*res.sqrtS);
        double p_d  = sqrt(std::max(E_d*E_d - m_d*m_d, 0.0));
        p_d_cm = Vec4D(E_d,  p_d*sinT*cos(phi),  p_d*sinT*sin(phi),  p_d*cosT);
        p_X_cm = Vec4D(E_pi,-p_d*sinT*cos(phi), -p_d*sinT*sin(phi), -p_d*cosT);
        if ( (res.part1->Flav().Kfcode()==2212 && res.part2->Flav().Kfcode()==2212 ) )
        {
            if ( res.isAnti ==0) fl_X   = Flavour(kf_pi_plus);
            else fl_X   = Flavour(kf_pi_plus).Bar();
        }
        else if ( res.part1->Flav().Kfcode()==2112 && res.part2->Flav().Kfcode()==2112 )
        {
            if (res.isAnti ==0) fl_X   = Flavour(kf_pi_plus).Bar();
            else fl_X   = Flavour(kf_pi_plus);
        }
        else fl_X   = Flavour(kf_pi);
    }
    else {
        // 3-body: d + pi + pi  — use Eq.(15), draw m_dpi2 uniformly
        double m_pi   = 0.139570;
        double m_dpi2_min = std::pow(m_d + m_pi, 2);
        double m_dpi2_max = std::pow(res.sqrtS - m_pi, 2);
        if (m_dpi2_max <= m_dpi2_min) { delete blob; return nullptr; }

        double m_dpi2 = m_dpi2_min + ran->Get()*(m_dpi2_max - m_dpi2_min);

        double E_d_term = (s + m_d*m_d - m_dpi2) / (2.0*res.sqrtS);  // the inner term
        double p_d      = sqrt(std::max(E_d_term*E_d_term - m_d*m_d, 0.0));  // Eq.(15)
        double E_d      = sqrt(p_d*p_d + m_d*m_d);  // recover E_d from p_d
        double E_pipi   = res.sqrtS - E_d;
        p_d_cm = Vec4D(E_d,   p_d*sinT*cos(phi),  p_d*sinT*sin(phi),  p_d*cosT);
        p_X_cm = Vec4D(E_pipi,-p_d*sinT*cos(phi), -p_d*sinT*sin(phi), -p_d*cosT);
        fl_X   = Flavour(kf_pi_plus);
    }

    // boost to lab
    boost.BoostBack(p_d_cm);
    boost.BoostBack(p_X_cm);

    Flavour fl_d(kf_deuterium);
    if (res.isAnti) fl_d = fl_d.Bar();

    Particle* deut   = new Particle(-1, fl_d, p_d_cm);
    Particle* recoil = new Particle(-1, fl_X,  p_X_cm);
    deut->SetPosition(res.part1->Position());
    recoil->SetPosition(res.part1->Position());
    blob->AddToOutParticles(deut);
    blob->AddToOutParticles(recoil);

    msg_Out() << "Deuteron formed via " << res.process
              << (res.isAnti ? " [ANTI]" : "")
              << " sqrtS=" << res.sqrtS << "\n";
    return blob;
}

Blob* Deuteron_Coalescence::operator()(Particle* A,
                                        Particle* B,
                                        const double& b2_fm2) {
    auto res = CheckPair(A, B);
    if (!res.part1) return nullptr;

    res = EvaluateFormation(res, b2_fm2);
    if (!res.formed) return nullptr;

    return BuildBlob(res);
}