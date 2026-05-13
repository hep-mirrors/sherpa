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

Blob* Deuteron_Coalescence::operator()(Particle* A,
                                        Particle* B,
                                        const double& b2_fm2) {
    Flavour flA = A->Flav(), flB = B->Flav();

    // identify p and n regardless of order
    Particle* proton  = nullptr;
    Particle* neutron = nullptr;
    bool isAnti = false;
    
    msg_Out()<<"Checking coalescence for particles: "<<flA<<" and "<<flB<<std::endl;

    if      (flA.Kfcode()==2212 && flB.Kfcode()==2112 &&
             !flA.IsAnti()      && !flB.IsAnti())
        { proton=A; neutron=B; isAnti=false; }
    else if (flA.Kfcode()==2112 && flB.Kfcode()==2212 &&
             !flA.IsAnti()      && !flB.IsAnti())
        { proton=B; neutron=A; isAnti=false; }
    else if (flA.Kfcode()==2212 && flB.Kfcode()==2112 &&
             flA.IsAnti()       && flB.IsAnti())
        { proton=A; neutron=B; isAnti=true;  }
    else if (flA.Kfcode()==2112 && flB.Kfcode()==2212 &&
             flA.IsAnti()       && flB.IsAnti())
        { proton=B; neutron=A; isAnti=true;  }
    else
    {
        msg_Out() << "Particles are not p+n or p~+n~, skipping coalescence"<<std::endl;
        return nullptr;
    }


    Vec4D p1 = proton->Momentum(); Vec4D p2 = neutron->Momentum();
    Vec4D P  = p1 + p2;
    double s     = P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3];
    double sqrtS = sqrt(std::max(s, 0.0));
    double m_d  = 1.875613;
    //m_deuteron.HadMass();
 

    //  relative momentum in CM frame 
    Poincare boost(P);
    Vec4D p1cm = p1;  boost.Boost(p1cm);
    double p_rel = sqrt(p1cm[1]*p1cm[1]+p1cm[2]*p1cm[2]+p1cm[3]*p1cm[3]);

    double Prob = 0.0;
    if (m_model == CoalescenceModel::HardDisk) {
        Prob = (b2_fm2 < m_b_coal*m_b_coal && p_rel < m_p_coal) ? 1.0 : 0.0;
    } 
    else {
        double C_b = (m_b_coal*m_b_coal) / (0.75*M_PI);
        double C_p = (m_p_coal*m_p_coal) / (0.75*M_PI);
        Prob = 0.75*std::exp(-b2_fm2/C_b) * std::exp(-(p_rel*p_rel)/C_p);
    }

    // this should be under an if, but i put here so i can see if this thing runs properly. 
    msg_Out() << "Coalescence candidate: "
              << (isAnti ? "p~+n~" : "p+n")
              << "  b="     << sqrt(b2_fm2) << " fm"
              << " b_coal =" << m_b_coal   << " fm"
              << " p_coal =" << m_p_coal   << " GeV"
              << "  p_rel=" << p_rel        << " GeV"
              << "  sqrt(s)=" << sqrtS      << " GeV"
              << "  m_d="   << m_d          << " GeV"
              << "  Prob="  << Prob         << "\n";

    if (Prob <= ran->Get()) return nullptr;

    if (sqrtS < m_d) {
        msg_Out() << "sqrt(s)=" << sqrtS << " < m_d=" << m_d
                  << ", cannot form deuteron\n";
        return nullptr;
    }
    
    Blob* blob = new Blob();
    blob->SetPosition(proton->Position());
    proton->SetStatus(part_status::decayed);
    neutron->SetStatus(part_status::decayed);
    blob->AddToInParticles(proton);
    blob->AddToInParticles(neutron);

    //  2-body kinematics: p+n → d+photon in CM frame 
    double E_d_cm = (s + m_d*m_d) / (2.0*sqrtS);
    double E_g_cm = (s - m_d*m_d) / (2.0*sqrtS);
    double p_d_cm = sqrt(std::max(E_d_cm*E_d_cm - m_d*m_d, 0.0));

    // isotropic photon direction in CM frame
    double cosT = 1.0 - 2.0*ran->Get();
    double sinT = sqrt(std::max(0.0, 1.0-cosT*cosT));
    double phi  = 2.0*M_PI*ran->Get();

    Vec4D p_d_cm_4(E_d_cm,
                    p_d_cm*sinT*cos(phi),
                    p_d_cm*sinT*sin(phi),
                    p_d_cm*cosT);
    Vec4D p_g_cm_4(E_g_cm,
                   -p_d_cm*sinT*cos(phi),
                   -p_d_cm*sinT*sin(phi),
                   -p_d_cm*cosT);

    // boost back to lab
    boost.BoostBack(p_d_cm_4);
    boost.BoostBack(p_g_cm_4);

    //  create outgoing particles 
    // Flavour fl_d   = isAnti ? m_deuteron.Bar() : m_deuteron;

    Flavour fl_d(kf_deuterium);

    Flavour fl_gam(kf_photon); 

    // msg_Out()<<"kf photon>"<< fl_gam.Kfcode() << std::endl;
    // msg_Out()<<"kf deuteron>"<< fl_d.Kfcode() << std::endl;

    //bottom lines problematic. 
    //   kf_deterium is problem.   
    Particle* deut  = new Particle(-1, fl_d,   p_d_cm_4);
    Particle* gamma = new Particle(-1, fl_gam,  p_g_cm_4);
    deut->SetPosition(proton->Position());
    gamma->SetPosition(proton->Position());

    blob->AddToOutParticles(deut);
    blob->AddToOutParticles(gamma);


    msg_Out() << "Deuteron formed: sqrt(s)=" << sqrtS
              << " E_gamma(CM)=" << E_g_cm << " GeV"
              << (isAnti ? " [ANTI]" : "") << "\n";

    return blob;
}