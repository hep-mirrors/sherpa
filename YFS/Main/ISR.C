#include "YFS/Main/ISR.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Model_Base.H"

#include <iostream>
#include <fstream>


#include <algorithm>
using namespace YFS;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;
// #define YFS_DEBUG


ISR::ISR()
{
  m_Kmin = sqrt(m_s) * (m_isrcut) / 2.;
  m_Kmax = sqrt(m_s * (1. - m_isrcut));
  m_cut = 1.; // set to 0 if anything fails to pass cuts
  m_nsuccess = 0;
  m_nfail = 0;
  m_ntotal = 0;
}

ISR::~ISR() {
}


double ISR::CalculateBeta(const Vec4D& p) {
  return (Vec3D(p).Abs() / p[0]);
}

void ISR::SetIncoming(YFS::Dipole *dipole) {
  p_dipole = dipole;
  m_beam1 = p_dipole->GetBornMomenta(0);
  m_beam2 = p_dipole->GetBornMomenta(1);
  m_b1 = CalculateBeta(m_beam1);
  m_b2 = CalculateBeta(m_beam2);
  m_mass = p_dipole->Mass();
  m_mass2 = m_mass * m_mass;
  m_am2 = sqr(m_beam1.Mass()+m_beam2.Mass()) / m_s;
  m_g  = p_dipole->m_gamma;
  m_gp = p_dipole->m_gammap;

  if (m_mass <= 0)THROW(fatal_error, "Initial leptons must be massive for YFS");

}

void ISR::NPhotons() {
  if (m_v < m_isrcut){
    m_n = 0;
    return;
  }
  m_nbar = m_gp * log(m_v / m_isrcut);
  if(m_fixed_photons!=-1) {
    m_n = m_fixed_photons;
  }
  if (m_nbar < 0 ) {
    msg_Error() << METHOD << "Warning: ISR photon average is less than 0" << std::endl;
  }
  int N = 0;
  double sum = 0;
  while (true) {
    N += 1;
    sum += log(ran->Get());
    if (sum <= -m_nbar) break;
  }
  m_n = N;
  if (m_n < 0) msg_Error() << METHOD << std::endl << "Nphotons < 0!!" << std::endl;
}



void ISR::GenerateAngles()
{
  // Generation of theta for two massive particles
  double weight = 1;
  if (m_kkmcAngles == 0) {
    double P = log((1.+m_b1)/(1.-m_b1))
                /(log((1.+m_b1)/(1.-m_b1))+log((1.+m_b2)/(1.-m_b2)));
    while (true) {
      if (ran->Get() < P) {
        double rnd = ran->Get();
        double a   = 1./m_b1*log((1.+m_b1)/(1.-m_b1));
        m_c        = 1./m_b1*(1.-(1.+m_b1)*exp(-a*m_b1*rnd));
        // m_c        = 1./m_b1*(pow(1+m_b1,rnd)/pow(1-m_b1,rnd-1)-1);
      }
      else {
        double rnd = ran->Get();
        double a   = 1./m_b2*log((1.+m_b2)/(1.-m_b2));
        m_c        = 1./m_b2*((1.-m_b2)*exp(a*m_b2*rnd)-1.);
        // m_c        = 1./m_b1*(pow(1-m_b2,rnd)/pow(1+m_b2,rnd-1)-1);
      }
      weight = 1.-((1.-m_b1*m_b1)/((1.-m_b1*m_c)*(1.-m_b1*m_c))
                        +(1.-m_b2*m_b2)/((1.+m_b2*m_c)*(1.+m_b2*m_c)))
                       /(2.*(1.+m_b1*m_b2)/((1.-m_b1*m_c)*(1.+m_b2*m_c)));
        if (ran->Get() < weight || m_kkmcAngles!=2) break;
      }
    if(m_kkmcAngles==2) m_angleWeight *= weight;
    m_theta = acos(m_c);
    m_sin = sin(m_theta);
    m_phi = 2.*M_PI * ran->Get();
    m_del1.push_back(1-m_b1*m_c);
    m_del2.push_back(1+m_b2*m_c);
    m_cos.push_back(m_c);
  }
  else {
    m_beta  = sqrt(1. - m_am2);
    double eps  = m_am2 / (1. + m_beta);
    double rn = ran->Get();
    double del1 = (2. - eps) * pow((eps / (2 - eps)), rn); // 1-beta*cos
    double del2 = 2. - del1;  // 1+beta*cos
    double costhg = (del2 - del1) / (2.*m_beta);
    double sinthg = sqrt(del1 * del2 - m_am2 * costhg * costhg);
    // symmetrization
    if (ran->Get() < 0.5) {
      double a = del1;
      del1 = del2;
      del2 = a;
      costhg = -costhg;
    }
    m_theta = acos(costhg);
    m_phi = 2.*M_PI * ran->Get();
    m_c = costhg;
    m_sin = sin(m_theta);
    m_del1.push_back(del1);
    m_del2.push_back(del2);
    m_cos.push_back(m_c);
  }
  if(abs(m_c)>1){
      msg_Error()<<"Photon angel out of bounds with cos(theta) = "<<m_c<<std::endl;
  }
}



void ISR::GeneratePhotonMomentum() {
  DEBUG_FUNC(METHOD << "# of soft photons to be generated: " << m_n << std::endl);
  Clean();
  if(m_v < m_isrcut && m_n!=0){
    msg_Error()<<"Warning: Generating real photon emissions below IR cut-off\n";
  }
  if(m_v > 1) msg_Error()<<"m_v > 1 = "<<m_v<<std::endl;
  if (m_n != 0) {
    Vec4D_Vector photons;
    GenerateAngles();
    double ms = m_mass2 / m_s;
    m_w = m_v;
    m_photon = {m_w,
                m_w * sin(m_theta) * cos(m_phi) ,
                m_w * sin(m_theta) * sin(m_phi) ,
                m_w * cos(m_theta)
               };
    m_photonSum += m_photon;
    m_photons.push_back(m_photon);

    double del1 = 1. - m_b1 * m_c;
    double del2 = 1. + m_b2 * m_c;
    // m_f    = m_alpi * (1 / ((del1) * (del2)) - ms / (del2 * del2) - ms / (del1 * del1));
    // m_fbar = m_alpi * (1. / ((del1) * (del2)));
    m_f = Eikonal(m_photon,m_beam1,m_beam2);
    m_fbar = EikonalMassless(m_photon,m_beam1,m_beam2);
    m_massW = m_f / m_fbar;
    m_yini.push_back(m_w * del1 / 2);
    m_zini.push_back(m_w * del2 / 2);
    for (int i = 1; i < m_n; i++) {
      GenerateAngles();
      m_w  = m_isrcut * pow(m_v / m_isrcut, ran->Get());
      del1 = 1. - m_b1 * m_c;
      del2 = 1. + m_b2 * m_c;
      // m_phi = 2.*M_PI * ran->Get();
      m_photon = { m_w,
                  m_w * sin(m_theta) * cos(m_phi) ,
                  m_w * sin(m_theta) * sin(m_phi) ,
                  m_w * cos(m_theta)
                 };
      m_photonSum += m_photon;
      m_photons.push_back(m_photon);
      // m_f    = m_alpi * (1. / ((del1) * (del2)) - ms / (del2 * del2) - ms / (del1 * del1));
      // m_fbar = m_alpi * (1. / ((del1) * (del2)));
      m_f = Eikonal(m_photon,m_beam1,m_beam2);
      m_fbar = EikonalMassless(m_photon,m_beam1,m_beam2);
      m_massW *= m_f / m_fbar;
      m_yini.push_back(m_w * del1 / 2);
      m_zini.push_back(m_w * del2 / 2);
    }
    MapPhotonMomentun();
  }
}


void ISR::MapPhotonMomentun() {
  if (m_n == 1) {
    m_diljac = 1.;
    m_lam = 1.;
    A = 0;
  }
  else {
    Vec4D P  = {2., 0., 0., 0. };
    double P2 = P * P;
    double K0 = m_photonSum[0];
    double PK = 2.*K0;
    double K2 = m_photonSum.Abs2();
    m_K2.push_back(K2);
    m_PTK.push_back(m_photonSum.PPerp());
    A = K2 * P2 /PK / PK;
    m_AA.push_back(A);
    m_lam  =  (m_v) * P2 / (PK) / (1. + sqrt(1. - A * m_v));
    if(IsBad(m_lam)){
      PRINT_VAR(m_v);
      PRINT_VAR(P2);
      PRINT_VAR(PK);
      PRINT_VAR(K2);
      PRINT_VAR(A);
      PRINT_VAR(m_photonSum);
      PRINT_VAR(m_b1);
      PRINT_VAR(m_b2);
    }
    m_lam0 = PK / P2 / m_v * (1. + sqrt(1. - m_v * A));
    m_scale.push_back(m_lam);
    m_diljac  = 0.5 * (1. + 1. / sqrt(1. - m_v * A));
  }
  
  m_diljac0 = 0.5 * (1. + pow(1. - m_v, -0.5));
  m_jacW = m_diljac / m_diljac0;
  m_jacvec.push_back(m_jacW);
  m_photonSum *= m_lam * sqrt(m_s) / 2.;

  for (size_t i = 0; i < m_photons.size(); ++i)
  {
    m_photons[i] *= m_lam * sqrt(m_s) / 2.;
    m_yini[i] /= m_lam;
    m_zini[i] /= m_lam;
    if (m_photons[i][0] < m_Kmin) m_cut = 0.;
  }
  if(m_photons.size()!=m_n){
    msg_Error()<<"Missmatch in Photon Multiplicity for ISR"<<std::endl
               <<" Poisson N = "<<m_n<<std::endl
               <<" Actual Photons = "<<m_photons.size()<<std::endl;
  }
}



void ISR::Weight() {
  m_ntotal += 1;
  double corrW = 1;
  if (m_v >= m_isrcut && m_n != 0 ) {
    m_weight = m_gp * pow(m_v, m_gp - 1) * m_diljac0 * pow(m_isrcut, m_g - m_gp);
  }
  else {
    m_massW = 1.0;
    m_jacW = 1.0;
    m_weight = m_g * pow(m_v, m_g - 1);
    double B = pow(m_isrcut, m_g) * (-m_g * m_isrcut + m_g + 1.) / (m_g + 1.);
    double D = pow(m_deltacut, m_g) * (-m_g * m_deltacut + m_g + 1.) / (m_g + 1.);
    corrW = 1. / (1. - D / B);
    m_weight *= corrW;
  }
  m_weight *= m_cut * m_massW * m_jacW * m_angleWeight;
  if (m_cut == 0) m_nfail += 1;
  else m_nsuccess += 1;
  DEBUG_FUNC("v = " << m_v << std::endl <<
             "vmin = " << m_isrcut << std::endl <<
             "vmax = " << m_vmax << std::endl <<
             "Kmin = " << m_Kmin << std::endl <<
             "Kmax = " << m_Kmax << std::endl <<
             "eps  = " << m_isrcut << std::endl <<
             "NPhotons  = " << m_n << std::endl <<
             "pow(m_v,m_gp-1) = " << pow(m_v, m_gp - 1) << std::endl <<
             "pow(m_isrcut,m_g-m_gp) = " << pow(m_isrcut, m_g - m_gp) << std::endl <<
             "Gamma = " << m_g << std::endl <<
             "Gamma prime = " << m_gp << std::endl <<
             "Lambda = " << m_lam << std::endl << "J0 = " << m_diljac0 << std::endl <<
             "J = " << m_diljac << std::endl << "W_J = " << m_jacW << std::endl << "cut = "
             << m_cut << std::endl
             << "Eikonal Weight = " << m_massW << std::endl
             << "Corrective Weight = " << corrW << std::endl
             << "Weight = " << m_weight << std::endl);
  for (size_t i = 0; i < m_photons.size(); ++i) msg_Debugging() << "k[" << i << "] = " << m_photons[i] << std::endl;
  if (IsBad(m_weight)) {
    msg_Error()<<METHOD<<std::endl << "YFS Weight is: " << m_weight << std::endl <<
                "v = " << m_v << std::endl <<
                "sqrt(m_sp) = " << sqrt(m_s*(1-m_v)) << std::endl <<
                "vmin = " << m_isrcut << std::endl <<
                "vmax = " << m_vmax << std::endl <<
                "weight = " << m_weight << std::endl <<
                "eps  = " << m_isrcut << std::endl <<
                "Nphotons  = " << m_n << std::endl <<
                "Gamma = " << m_g << std::endl <<
                "Gamma prime = " << m_gp << std::endl <<
                "Kmin = " << m_Kmin << std::endl <<
                "Lambda = " << m_lam << std::endl << "J0 = " << m_diljac0 << std::endl <<
                "J = " << m_diljac << std::endl << "W_J = " << m_jacW << std::endl << "cut = "
                << m_cut << std::endl
                << "W_mass = " << m_massW << std::endl;
  }
}


void ISR::Clean() {
  m_yini.clear();
  m_zini.clear();
  m_del1.clear();
  m_del2.clear();
  m_photons.clear();
  m_photonSum = Vec4D(0, 0, 0, 0);
  m_weight = 1.0;
  m_angleWeight = 1.0;
  m_massW = 1.0;
  m_jacW  = 1.0;
  m_cut   = 1.0;
  m_lam0 = 1.0;
  m_diljac0 = 1.0;
  m_diljac = 1.0;
  m_cos.clear();
  m_scale.clear();
  m_jacvec.clear();
  m_AA.clear();
  m_K2.clear();
  m_PTK.clear();
}



void ISR::MakeYFS() {
  GeneratePhotonMomentum();
}

double ISR::Flux() {
  double m_f = (m_beam1 + m_beam2).Abs2() / (m_qvec1 + m_qvec2).Abs2();
  return m_f;
}

void ISR::Sort(Vec4D_Vector & p) {
  struct {
    bool operator()(Vec4D a, Vec4D b) const
    {
      return a.PPerp() > b.PPerp();
    }
  } cmp;

  std::sort(p.begin(), p.end(), cmp);
}

void ISR::Sort(std::vector<double> &p) {
  struct {
    bool operator()(double a, double b) const
    {
      return a < b;
    }
  } cmp;

  std::sort(p.begin(), p.end(), cmp);
}

double ISR::Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}

double ISR::EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  // return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
  return m_alpha/(4*M_PI*M_PI)*(2*p1*p2/((p1*k)*(p2*k)));
}