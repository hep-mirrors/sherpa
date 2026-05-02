#include "HADRON_RESCATTERING/Main/Hadron_Rescatterings.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescatterings::Hadron_Rescatterings(const bool & on) :
  m_on(on), m_nfails(0), m_mm2mb(1.e25)
{
  hrpars = new HR_Parameters();
}

void Hadron_Rescatterings::Initialize() {
  // ScatteringBase* p_scatteringBase = NULL;
  m_mm = new MesonMeson();
  m_bm = new BaryonMeson();
  m_bb = new BaryonBaryon();
}

Hadron_Rescatterings::~Hadron_Rescatterings() {
  delete m_mm;
  delete m_bm;
  delete m_bb;
  if (hrpars) { delete hrpars; hrpars = NULL; }
}

Blob * Hadron_Rescatterings::operator()(Particle * A,Particle * B,
					const double & dist2) {

  Vec4D p1 = A->Momentum();
  Vec4D p2 = B->Momentum();
  Vec4D P  = p1 + p2;
  double comEnergy = P[0]*P[0]- P[1]*P[1]-P[2]*P[2]-P[3]*P[3];
  ScatteringBase* interaction = nullptr;
  const Flavour flA = A->Flav(), flB = B->Flav();

  if (flA.IsMeson() && flB.IsMeson())         interaction = m_mm;
  else if (flA.IsBaryon() && flB.IsBaryon())  interaction = m_bb;
  else                                        interaction = m_bm;
  
  double xstot = interaction->XStot(flA,flB,comEnergy); // return value is in millibarns.
  xstot = 0.1*xstot; // millibarn to fm^2. 
  double C = xstot/(0.75 * M_PI);
  double dist2_fm2 = dist2 * 1e24;          // mm^2 → fm^2
  double r = ATOOLS::ran->Get();
  double Prob = std::exp(-dist2_fm2/C);
  msg_Out() << "dist2(mm^2)=" << dist2 
          << " dist2(fm^2)=" << dist2_fm2
          << " dist(fm)="   << sqrt(dist2_fm2)
          << " C(fm^2)="     << C << std::endl;

  msg_Out()<<"Prob is : "<< Prob << std::endl; 
  msg_Out()<<"random value is: "<< r<<std::endl;
  if (Prob > r){
    msg_Out()<<"Rescatter maybe: "
    <<METHOD<<"("<<flA<<", "<<flB<<"):\n; "
    <<"xstot = "<<xstot<<" vs. dist^2 = "<<dist2 
    <<" --> "<<(M_PI*dist2*m_mm2mb)<<"Prob = "<<Prob<<std::endl;

    Poincare boost(P);
    Vec4D p1cm = A->Momentum();  boost.Boost(p1cm);   // modifies in place
    Vec4D p2cm = B->Momentum();  boost.Boost(p2cm);
    msg_Out()<< "p1 in CM frame: " << p1cm[1]<<", "<<p1cm[2]<<", "<<p1cm[3]<<std::endl;
    msg_Out()<< "p2 in CM frame: " << p2cm[1]<<", "<<p2cm[2]<<", "<<p2cm[3]<<std::endl;
    //after boost along with P, the total momentum, they are back to back. 

    // magnitude of 3-momentum in CM frame. this should be equal to p_2cm^2, since after the boost, they are back to back.
    double p_cm = sqrt(p1cm[1]*p1cm[1] + p1cm[2]*p1cm[2] + p1cm[3]*p1cm[3]);
    //Elastic scattering:: m_1 == m_3, m_2 == m_4, but with new directions.
    // then momentum transfer t = (p1-p3)^2 = (p1cm-p3cm)^2, since E_1 == E_3.
    // Then t = -2 p_cm^2 (1-cos(theta)), where theta is the scattering angle in CM frame.
    // Assume an exponential distribution in t, i.e. P ~ exp(b t) = exp(-2 b p_cm^2 (1-cos(theta))), where b is the slope parameter.
    // get cosθ from e^{b|t|}, t = -2p_cm²(1-costheta), since P ~ e^{b t} = e^{-2 b p_cm²(1-costheta)}. 
    
    double bSlope = 2.4; 
    //why is bSlope = 2.4 ? https://arxiv.org/pdf/hep-ph/0312187 --> page 5, eq 4a.
    double lam     = 2.0 * bSlope * p_cm * p_cm;
    double r2      = ran->Get();
    double cos_theta; 
    //consider small and large angle. 
    if (lam < 1e-8)
        cos_theta = 1.0 - 2.0 * r2;              
    else
        cos_theta = 1.0 + log(1.0 - r2*(1.0 - exp(-2.0*lam))) / lam;

    double sin_theta = sqrt(std::max(0.0, 1.0 - cos_theta*cos_theta));
    double phi       = 2.0 * M_PI * ran->Get();

    //  building q1,q2 in CM frame
    // Same energy as incoming (elastic), but with a new direction. angles are like this so that the magnitudes are eq. 
    double qx = p_cm * sin_theta * cos(phi);
    double qy = p_cm * sin_theta * sin(phi);
    double qz = p_cm * cos_theta;

    Vec4D q1cm(p1cm[0],  qx,  qy,  qz);
    Vec4D q2cm(p2cm[0], -qx, -qy, -qz);   // back-to-back in CM frame again. 

    Vec4D dp = p1cm - q1cm;
    double t = dp[0]*dp[0] - dp[1]*dp[1] - dp[2]*dp[2] - dp[3]*dp[3];
    msg_Out() << "Mandelstam t = " << t << "  (should be <= 0)"<<std::endl;

    // boost q1,q2 back to lab frame 
    // BoostBack() is the exact inverse of Boost()
    Vec4D q1lab = q1cm;  boost.BoostBack(q1lab);
    Vec4D q2lab = q2cm;  boost.BoostBack(q2lab);

    // ── sanity check 
    msg_Out() << "Energy conservation: "
              << (q1lab[0]+q2lab[0]) << " vs " << (p1[0]+p2[0]) << std::endl;

    Blob* rblob = new Blob();
    rblob->SetType(btp::Elastic_Collision);
    rblob->SetPosition(A->Position());

    //they gotta decay, but not sure if status matters at the moment. 
    A->SetStatus(part_status::decayed);
    B->SetStatus(part_status::decayed);
    rblob->AddToInParticles(A);
    rblob->AddToInParticles(B);

    Particle* pout1 = new Particle(-1, flA, q1lab);
    Particle* pout2 = new Particle(-1, flB, q2lab);
    pout1->SetPosition(pout1->Position());
    pout2->SetPosition(pout2->Position());
    rblob->AddToOutParticles(pout1);
    rblob->AddToOutParticles(pout2);     
    return rblob; 
  }

  
    return NULL; 
}

void Hadron_Rescatterings::Reset() {}

