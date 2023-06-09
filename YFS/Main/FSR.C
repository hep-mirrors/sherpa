
#include "YFS/Main/FSR.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include <iostream>
#include <fstream>


#include <algorithm>
using namespace YFS;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;
using namespace PHASIC;

ofstream myfile;

FSR::FSR()
{
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  s["FSR_EMIN"].SetDefault(1e-9);
  s["FSR_FCUT"].SetDefault(0);
  s["FSR_NBAR"].SetDefault(0);
  s["MASSIVE_NBAR"].SetDefault(0);
  s["FSR_EIK"].SetDefault(0);
  s["FSR_CRU"].SetDefault(1);
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  s["FSR_CUT"].SetDefault(1e-3);
  m_Edelta = s["FSR_EMIN"].Get<double>();
  m_kkmcAngles = s["KKMC_ANG"].Get<bool>();
  m_fsrcut = s["FSR_CUT"].Get<double>() * m_vmin;
  m_fsrcutF = s["FSR_FCUT"].Get<double>();
  m_nbar = s["FSR_NBAR"].Get<double>();
  m_use_massive_nbar = s["MASSIVE_NBAR"].Get<bool>();
  m_use_crude = s["FSR_CRU"].Get<int>();
  m_eikonal_mode = s["FSR_EIK"].Get<int>();
  p_fsrFormFact = new YFS::YFS_Form_Factor();
  m_omegaMin = m_fsrcut;
  m_omegaMax = std::numeric_limits<double>::max();
  // m_alpi = m_alpha / M_PI;
}

FSR::~FSR() {
  if(p_yfsFormFact) delete p_yfsFormFact;
  delete p_fsrFormFact;
}

bool FSR::Initialize(YFS::Dipole_Vector::iterator dipole) {
  p_dipole = &dipole[0];
  m_fsrWeight = 0.;
  m_mass.clear();
  m_dipole.clear();
  m_dipoleFl.clear();
  m_dipole.push_back(p_dipole->GetOldMomenta(0));
  m_dipole.push_back(p_dipole->GetOldMomenta(1));
  m_dipoleFl.push_back(p_dipole->GetFlav(0));
  m_dipoleFl.push_back(p_dipole->GetFlav(1));
  m_QFrame = m_dipole[0] + m_dipole[1];
  for (size_t i = 0; i < m_dipole.size(); ++i) m_mass.push_back(p_dipole->m_masses[i]);
  if (IsZero(m_mass[0]) || IsZero(m_mass[1])) {
    THROW(fatal_error, "Charged particles must be massive for YFS");
  }
  m_Q1 = p_dipole->m_charges[0];
  m_Q2 = p_dipole->m_charges[1];
  m_QF2 = m_Q1 * m_Q2;
  m_dip_sp = p_dipole->Sprime();
  MakePair(sqrt(m_dip_sp), m_bornQ1, m_bornQ2);
  m_EQ = sqrt(m_dip_sp) / 2.;
  m_Emin = 0.5 * sqrt(m_s) * m_vmin;
  m_Kmin = 0.5 * m_fsrcut * sqrt(m_dip_sp);
  m_Kmax = sqrt(m_dip_sp) / 2.;
  m_hideW = 1.;
  if (m_dipole.size() != 2) {
    THROW(fatal_error, "Dipole size incorrect in YFS FSR")
  }
  m_p1p2 = m_dipole[0] * m_dipole[1];
  m_beta1 = CalculateBeta(m_dipole[0]);
  m_beta2 = CalculateBeta(m_dipole[1]);
  m_mu1 = 1. - sqr(m_beta1);
  m_mu2 = 1. - sqr(m_beta2);
  m_g  = p_dipole->m_gamma;
  m_gp = p_dipole->m_gammap;
  if (m_use_massive_nbar) m_nbar = -m_g * log(m_fsrcut);
  else m_nbar = -m_gp * log(m_fsrcut);

  if (IsBad(m_nbar)) {
    PRINT_VAR(m_dipole);
    PRINT_VAR(m_g);
    PRINT_VAR(m_gp);
    PRINT_VAR(m_mass);
    PRINT_VAR(m_QF2);
    PRINT_VAR(m_betaf);
    PRINT_VAR(m_amc2);
    PRINT_VAR(m_dip_sp);
  }
  p_fsrFormFact->SetCharge((-m_QF2));
  return true;
}




double FSR::CalculateBeta(const Vec4D& p) {
  return Vec3D(p).Abs() / p[0];
}

void FSR::CalculateBetaBar() {
  m_betaBar = m_dip_sp - sqr(m_mass[0] - m_mass[1]);
  m_betaBar *= m_dip_sp - sqr(m_mass[0] + m_mass[1]);
  m_betaBar = sqrt(m_betaBar) / m_dip_sp;
  m_betaBar1 = CalculateBeta(m_dipole[0]);
  m_betaBar2 = CalculateBeta(m_dipole[1]);
  m_betaBar1 = m_betaBar2 = m_betaBar;
  double t1 = (1. + m_betaBar1 * m_betaBar2) / (m_betaBar1 + m_betaBar2);
  double logarg =  2.*(1. + m_betaBar1) * (1. + m_betaBar2) / ((1. - m_betaBar1) * (1. - m_betaBar2));
  m_gBar  = -m_QF2 * m_alpi * t1 * (log(logarg / 2.) - 2.); // See Mareks phd thesis A.2.1
  m_gpBar = -m_QF2 * m_alpi * t1 * (log(logarg / 2.));

  if (IsBad(m_betaBar)) {
    PRINT_VAR(m_sX);
    PRINT_VAR(sqr(m_mass[0] + m_mass[1]));
    PRINT_VAR(m_mass);
  }
}

void FSR::GenerateAngles() {
// Generation of theta for two massive particles
  double del1, del2;
  double am2 = 4.*sqr((m_mass[0])) / m_dip_sp;

  double P  = log((1. + m_beta1) / (1. - m_beta1))
              / (log((1. + m_beta1) / (1. - m_beta1)) + log((1. + m_beta2) / (1. - m_beta2)));
  if (!m_kkmcAngles) {
    while (true) {
      m_c = 0.;
      if (ran->Get() < P) {
        double rnd = ran->Get();
        double a   = log((1. + m_beta1) / (1. - m_beta1));
        m_c        = 1. / m_beta1 * (1. - (1. + m_beta1) * exp(-a * rnd));
      }
      else {
        double rnd = ran->Get();
        double a   = log((1. - m_beta2) / (1. + m_beta2));
        m_c        = -1. / m_beta2 * (1. - (1. - m_beta2) * exp(-a * rnd));
      }
      double weight = 1. - ((1. - m_beta1 * m_beta1) / sqr(1. - m_beta1 * m_c)
                            + (1. - m_beta2 * m_beta2) / sqr(1. + m_beta2 * m_c))
                      / (2.*(1. + m_beta1 * m_beta2) / ((1. - m_beta1 * m_c) * (1. + m_beta2 * m_c)));
      if (ran->Get() < weight) break;
      if (weight < 0) {
        msg_Error() << "FSR angle acceptance weight is < 0.\n Accepting Cos(theta) as " << m_c << std::endl;
        break;
      }
    }
    m_theta = acos(m_c);
    m_phi   = 2.*M_PI * ran->Get();
    m_st = 1. - sqr(m_c);
    del1 = 1 - m_beta1 * m_c;
    del2 = 1 + m_beta2 * m_c;
    // PRINT_VAR(m_theta);
  }
  else {
    double beta  = sqrt(1. - am2);
    double eps  = am2 / (1. + beta);
    double rn = ran->Get();                    // 1-beta
    del1 = (2. - eps) * pow((eps / (2 - eps)), rn); // 1-beta*costhg
    del2 = 2. - del1;  // 1+beta*costhg
    // calculation of sin and cos theta from internal variables
    double costhg = (del2 - del1) / (2.*beta);         // exact
    double sinthg = sqrt(del1 * del2 - am2 * costhg * costhg); // exact
    // symmetrization
    if (ran->Get() < 0.5) {
      double a = del1;
      del1 = del2;
      del2 = a;
      costhg = -costhg;
    }
    // double weight = 1. - ((1. - m_beta1 * m_beta1) / sqr(1. - m_beta1 * m_c)
    //                       + (1. - m_beta2 * m_beta2) / sqr(1. + m_beta2 * m_c))
    //                 / (2.*(1. + m_beta1 * m_beta2) / ((1. - m_beta1 * m_c) * (1. + m_beta2 * m_c)));
    m_theta = acos(costhg);
    m_phi = 2.*M_PI * ran->Get();
    m_c = costhg;
    m_st = sinthg;
    del1 = 1 - m_beta1 * m_c;
    del2 = 1 + m_beta2 * m_c;
  }
  m_cos.push_back(m_c);
  m_sin.push_back(m_st);
  m_fbarvec.push_back(1. / (del1 * del2) * (1. - am2 / 2.));
}



void FSR::NPhotons() {
  int N = 0;
  double sum = 0.0;
  while (true) {
    N += 1;
    sum += log(ran->Get());
    m_rand.push_back(sum / (-m_nbar));
    if (sum <= -m_nbar) break;
  }
  m_n = N - 1;
  m_N = m_n;
  p_dipole->SetNPhoton(m_N);
  if (m_n < 0) msg_Error() << METHOD << std::endl << "Nphotons < 0!!" << std::endl;
}


void FSR::GeneratePhotonMomentum() {
  // DEBUG_FUNC(METHOD<<"# of soft photons generated: "<<m_n<<std::endl);
  DEBUG_FUNC(" FSR Nphotons: " << m_n);
  m_photons.clear();
  m_MassWls.clear();
  m_massW = 1.0;
  m_photonSum = Vec4D(0, 0, 0, 0);
  m_cos.clear();
  m_sin.clear();
  m_yini.clear();
  m_zini.clear();
  m_k0.clear();
  m_fbarvec.clear();
  m_rand.clear();
  m_phi_vec.clear();
  for (int i = 0; i < m_n; i++) {
    GenerateAngles();
    double k0 = pow(m_fsrcut, ran->Get());
    Vec4D photon = {k0,
                    k0 * m_st * cos(m_phi) ,
                    k0 * m_st * sin(m_phi) ,
                    k0 * m_c
                   };
    m_photons.push_back(photon);
    m_photonSum += photon;
    m_k0.push_back(k0);
  }
}


bool FSR::MakeFSR() {
  m_dist1.clear();
  m_dist2.clear();
  m_del1.clear();
  m_del2.clear();
  m_photonspreboost.clear();
  m_cut = 1;
  m_wt2 = 1.0;
  m_yy = 1.0;
  m_xfact = 1.0;
  m_sQ = m_dip_sp;
  double smin = (m_mass[0] + m_mass[1]);
  NPhotons();
  GeneratePhotonMomentum(); // Run this even if no photons to clear previous event
  if (m_photons.size() == 0) {
    m_sprim = m_dip_sp;
    m_sQ = m_sprim;
  }
  else {
    if (m_photonSum.E() >= 1) {
      RejectEvent();
      m_cut = 2;
      return false;
    }
    RescalePhotons();
    m_sQ = m_dip_sp * m_yy;
    if ( sqrt(m_sQ) < smin ) {
      RejectEvent();
      m_cut = 3;
      return false;
    }
  }
  if (m_n == 0 && m_wt2 != 1) {
    msg_Error() << METHOD << "Incorrect jacobian in YFS FSR" << std::endl;
  }
  m_u = 1 - m_sprim / m_dip_sp;
  MakePair(sqrt(m_sprim), m_dipole[0], m_dipole[1]);
  m_px = m_dipole[0] + m_dipole[1] + m_photonSum;
  m_sX = m_px.Abs2();
  m_Q = m_dipole[0] + m_dipole[1];
  double masc1 = m_mass[0] * sqrt(m_sQ / m_dip_sp);
  double masc2 = m_mass[1] * sqrt(m_sQ / m_dip_sp);
  CE.Isotropic2Momenta(m_Q, sqr(masc1), sqr(masc2), m_r1, m_r2, ran->Get(), ran->Get(), -1, -1);
  CalculateBetaBar();
  p_dipole->AddToGhosts(m_r1);
  p_dipole->AddToGhosts(m_r2);
  for (int i = 0; i < 2; ++i) p_dipole->SetMomentum(i, m_dipole[i]);
  m_photonSumPreBoost = m_photonSum;
  if (m_cut != 0) return true;
  else return false;
}

void FSR::RescalePhotons() {
  m_xfact = 1. / (1. - m_photonSum[0]);
  for (int i = 0; i < m_photons.size(); ++i) m_photons[i] *= m_xfact;
  m_photonSum *= m_xfact;
  // K = m_photonSum;
  m_yy = 1. / (1. + m_photonSum[0] + 0.25 * m_photonSum.Abs2());
  m_wt2 = m_yy * (1. + m_photonSum[0]);
  m_sprim = m_dip_sp * m_yy;
  m_sqtest = m_sprim;

  // Rescale all photons
  double ener = sqrt(m_sprim) / 2.;
  for (size_t i = 0; i < m_photons.size(); ++i) {
    m_photons[i] *= ener;
    m_photonspreboost.push_back(m_photons[i]);
  }
  p_dipole->AddPhotonsToDipole(m_photons);
  m_photonSum *= ener;
  for (auto k : m_photons) {
    msg_Debugging() << k << std::endl;
  }
}

bool FSR::F(ATOOLS::Vec4D_Vector &k) {
  double del1, del2;
  m_sqtest = m_dip_sp;
  m_sX = (m_dipole[0] + m_dipole[1] + m_photonSum).Abs2();
  double ener = sqrt(m_sprim) / 2.;
  double am1 = sqr(m_dipole[0].Mass() / ener);
  double am2 = sqr(m_dipole[1].Mass() / ener);
  m_eta1 = (m_sprim + m_dipole[0].Abs2() - m_dipole[1].Abs2()) / m_sprim;
  m_eta2 = (m_sprim - m_dipole[0].Abs2() + m_dipole[1].Abs2()) / m_sprim;
  Vec4D p1 = m_dipole[0];
  Vec4D p2 = m_dipole[1];
  double eta1(0), eta2(0);
  double betan = sqrt((m_sprim - pow(m_mass[0] - m_mass[1], 2)) * (m_sprim - pow(m_mass[0] + m_mass[1], 2))) / m_sprim;
  MakePair(sqrt(m_sprim), p1, p2, p1.Mass(), p2.Mass() , eta1, eta2);
  CalculateBetaBar();
  for (size_t i = 0; i < k.size(); ++i)
  {
    if (m_cos[i] <= 0.) {
      del1 = am1 / (m_eta1 + betan) + betan * sqr(m_sin[i]) / (1. + m_cos[i]);
      del2 = m_eta2 + betan * m_cos[i];
    }
    else {
      del1 = m_eta1 - betan * m_cos[i];
      del2 = am2 / (m_eta2 + betan) + betan * sqr(m_sin[i]) / (1. - m_cos[i]);
    }
    m_del1.push_back(del1);
    m_del2.push_back(del2);
    if (m_eikonal_mode == 1) {
      m_f = Eikonal(k[i]);
      m_fbar = EikonalInterferance(k[i]);
    }
    else {
      m_f = 1. - (am1 + am2) / 4 - am1 / 4.*del2 / del1 - am2 / 4.*del1 / del2;
      m_f /= del1 * del2;
      m_fbar = m_fbarvec[i];
    }
    if (IsBad(m_f)) {
      PRINT_VAR(m_fbar);
      PRINT_VAR(del1);
      PRINT_VAR(del2);
      PRINT_VAR(betan);
      PRINT_VAR(sqrt(m_sQ));
      m_f = 0;
    }
    m_MassWls.push_back(m_f / m_fbar);
    m_dist1.push_back(m_f);
    m_dist2.push_back(m_fbar);
    if (IsBad(m_massW)) {
      PRINT_VAR(m_f);
      PRINT_VAR(m_fbar);
      PRINT_VAR(m_betaBar);
      PRINT_VAR(m_cos[i]);
      return false;
    }
  }
  return true;
}

void FSR::YFS_FORM() {
  // r1, r2 are the corresponding q* vectors defined pg 46 arxiv  9912214
  // they are created such that sqr(r_i) = sqr(m_i)sQ/sX
  // create back to back
  for (int i = 0; i < 2; ++i) m_dipole[i] = p_dipole->GetMomenta(i);
  m_photons = p_dipole->GetPhotons();
  m_photonSum = p_dipole->GetPhotonSum();
  m_r1 = p_dipole->GetGhost(0);
  m_r2 = p_dipole->GetGhost(1);
  m_Q = m_dipole[0] + m_dipole[1];
  m_sQ = m_Q.Abs2();
  // PRINT_VAR(sqrt(m_sX));
  m_sX = (m_dipole[0] + m_dipole[1] + m_photonSum).Abs2();
  CalculateBetaBar();
  m_Q = m_dipole[0] + m_dipole[1];
  m_q1q2 = m_dipole[0] * m_dipole[1];
  double Eqq = 0.5 * sqrt(m_sQ);
  double Delta1 = m_fsrcut * (1. + 2.*m_Q * m_photonSum / m_sQ);
  m_delta1 = Delta1;
  m_EminQ  = Eqq * m_fsrcut * (1. + 2.*m_Q * m_photonSum / m_sQ);
  m_expf = exp(m_volmc);
  double mass1sc = m_mass[0] * sqrt(m_sQ / m_sX);
  double mass2sc = m_mass[1] * sqrt(m_sQ / m_sX);
  double Eq1   = (m_sQ + m_mass[0] * m_mass[0] - m_mass[1] * m_mass[1]) / (2 * sqrt(m_sQ));
  double Eq2   = (m_sQ + m_mass[1] * m_mass[1] - m_mass[0] * m_mass[0]) / (2 * sqrt(m_sQ));
  double EQQ = sqrt(m_sQ) / 2.;
  m_bvrA = p_fsrFormFact->A(m_q1q2, m_mass[0], m_mass[1]);
  double YFS_IR = -2.*m_alpi * abs(m_QF2) * (m_q1q2 * p_fsrFormFact->A(m_q1q2, m_mass[0], m_mass[1]) - 1.) * log(1 / m_delta1);

  m_btilStar = p_fsrFormFact->BVR_full(m_q1q2, m_dipole[0][0], m_dipole[1][0], m_mass[0], m_mass[1], m_Emin, m_photonMass, 0);
  m_btil     = p_fsrFormFact->BVR_full(m_q1q2, Eq1, Eq2, m_mass[0], m_mass[1], m_EminQ, m_photonMass, 0);

  if (m_use_crude) {
    m_BtiXcru = p_fsrFormFact->BVR_cru(m_r1 * m_r2, m_r1[0], m_r2[0], m_r1.Mass(), m_r2.Mass(), m_Emin, m_photonMass);
    m_BtiQcru = p_fsrFormFact->BVR_cru(m_r1 * m_r2, Eqq, Eqq, m_r1.Mass(), m_r2.Mass(), m_EminQ, m_photonMass);
  }
  else {
    m_BtiXcru = p_fsrFormFact->BVR_full(m_r1 * m_r2, m_r1[0], m_r2[0], m_r1.Mass(), m_r2.Mass(), m_Emin, m_photonMass, 0);
    m_BtiQcru = p_fsrFormFact->BVR_full(m_r1 * m_r2, Eqq, Eqq, m_r1.Mass(), m_r2.Mass(), m_EminQ, m_photonMass, 0);
  }
  m_A4 = p_fsrFormFact->A4(m_r1 * m_r2, m_r1.E(), m_r2.E(), m_r1.Mass(), m_r2.Mass());
  m_A = p_fsrFormFact->A(m_r1 * m_r2, m_r1.Mass(), m_r2.Mass());
  m_DelYFS = m_btilStar - m_btil;
  m_delvol = m_BtiXcru  - m_BtiQcru;
  double amc = 4 * sqr(m_mass[0]) / m_sX;
  double bb = sqrt(1 - amc);
  double L = 2 * m_alpha / M_PI * (1 + bb * bb) / (2 * bb) * log(sqr(1 + bb) / amc);
  m_volmc = m_gpBar * log(1. / m_fsrcut) - m_delvol;
  m_hideW = exp(YFS_IR + m_DelYFS + m_volmc);
  m_YFS_IR = exp(YFS_IR + m_DelYFS);
}

void FSR::HidePhotons() {
  m_NRemoved  = 0;
  m_massW = 1.;
  m_yini.clear();
  m_zini.clear();
  Vec4D_Vector ph;
  std::vector<int> mark;
  std::vector<double> y, z, del1, del2;
  m_photonSum *= 0;
  m_photons = p_dipole->GetPhotons();
  // mark photons for removal
  for (size_t i = 0; i < m_photons.size(); ++i) {
    if (m_photons[i].E() <= m_Emin) {
      m_NRemoved++;
      mark.push_back(i);
    }
    else {
      m_massW *= m_MassWls[i];
      ph.push_back(m_photons[i]);
      del1.push_back(m_del1[i]);
      del2.push_back(m_del2[i]);
    }
  }

  for (size_t i = 0; i < ph.size(); ++i)
  {
    m_photonSum += ph[i];
    m_yini.push_back(ph[i].E()*del1[i] / sqrt(m_sprim));
    m_zini.push_back(ph[i].E()*del2[i] / sqrt(m_sprim));

    //m_photons.erase(m_photons.begin()+mark[i]);
  }
  p_dipole->SetNPhoton(ph.size());
  m_photons = ph;
  p_dipole->AddPhotonsToDipole(m_photons);
  
}



void FSR::MakePair(double cms, Vec4D &p1, Vec4D &p2, double mass1, double mass2,
                   double &eta1, double &eta2) {
  double E = cms / 2.;
  double s = sqr(cms);
  double beta2 = (s - sqr(mass1 - mass2)) * (s - sqr(mass1 + mass2)) / (s * s);
  double beta =  sqrt(beta2);
  eta1 = (s + sqr(mass1) - sqr(mass2)) / (s);
  eta2 = (s - sqr(mass1) + sqr(mass2)) / (s);
  p1 = {eta1 * E, 0, 0, beta * E};
  p2 = {eta2 * E, 0, 0, -beta * E};
  if (!IsEqual(p1.Mass(), mass1, 1e-5) || !IsEqual(p2.Mass(), mass2, 1e-5)) {
    msg_Error() << METHOD << "Error in masses for energy = " << cms << std::endl
                << "s = " << s << std::endl
                << "beta2 = " << beta2 << std::endl
                << "beta = " << beta << std::endl
                << "E = " << E << std::endl
                << "Mass of p1 = " << p1.Mass() << std::endl
                << "p1 = " << p1 << std::endl
                << "Mass should be = " << mass1 << std::endl
                << "Difference = " << abs(p1.Mass() - mass1) << std::endl
                << "Mass of p2 = " << p2.Mass() << std::endl
                << "p2 = " << p2 << std::endl
                << "Mass should be = " << mass2 << std::endl
                << "Difference = " << abs(p2.Mass() - mass2) << std::endl;
  }
}

void FSR::MakePair(double cms, Vec4D &p1, Vec4D &p2) {
  double E = cms / 2.;
  double s = sqr(cms);
  double mass1 = p1.Mass();
  double mass2 = p2.Mass();
  double beta2 = (s - sqr(mass1 - mass2)) * (s - sqr(mass1 + mass2)) / (s * s);
  double beta =  sqrt(beta2);
  double eta1 = (s + sqr(mass1) - sqr(mass2)) / s;
  double eta2 = (s - sqr(mass1) + sqr(mass2)) / s;
  p1 = {E * eta1, 0, 0, beta * E};
  p2 = {E * eta2, 0, 0, -beta * E};
  if (!IsEqual(p1.Mass(), mass1, 1e-5) || !IsEqual(p2.Mass(), mass2, 1e-5)) {
    msg_Error() << METHOD << "Error in masses for energy = " << cms << std::endl
                << "s = " << s << std::endl
                << "beta2 = " << beta2 << std::endl
                << "beta = " << beta << std::endl
                << "E = " << E << std::endl
                << "Mass of p1 = " << p1.Mass() << std::endl
                << "p1 = " << p1 << std::endl
                << "Mass should be = " << mass1 << std::endl
                << "Difference = " << abs(p1.Mass() - mass1) << std::endl
                << "Mass of p2 = " << p2.Mass() << std::endl
                << "p2 = " << p2 << std::endl
                << "Mass should be = " << mass2 << std::endl
                << "Difference = " << abs(p2.Mass() - mass2) << std::endl;
  }
}

void FSR::BoostDipole(Vec4D_Vector &dipole) {
  Vec4D QMS = dipole[0] + dipole[1] + m_photonSum;
  ATOOLS::Poincare poin(QMS);
  poin.Boost(dipole[0]);
  poin.Rotate(dipole[0]);
  poin.Boost(dipole[1]);
  poin.Rotate(dipole[1]);
}

void FSR::Weight() {
  CalculateBetaBar();
  // m_fsrWeight=m_massW*m_hideW*m_wt2;
  double spp = m_sprim * m_sX / m_dip_sp;
  double form = p_fsrFormFact->BVV_full(m_dipole[0], m_dipole[1], 1e-15, sqrt(m_sQ) / 2., 1);
  if (IsBad(form)) form = exp(m_gBar / 4 + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5));
  if (m_photons.size() == 0) m_sprim = m_sX;

  double L = 2 * m_alpha / M_PI * (log(m_sprim / sqr(m_mass[0]) - 1));
  if (sqrt(m_s) > 150) form = m_gBar;
  else form = L;
  // PRINT_VAR(m_dip_sp/m_sQ);
  // PRINT_VAR(m_mass);
  m_betaBar = m_sprim - sqr(m_mass[0] - m_mass[1]);
  m_betaBar *= m_sprim - sqr(m_mass[0] + m_mass[1]);
  // double amc = sqr((m))
  m_betaBar = sqrt(m_betaBar) / m_sprim;
  m_betaBar1 = CalculateBeta(m_dipole[0]);
  m_betaBar2 = CalculateBeta(m_dipole[1]);
  m_betaBar1 = m_betaBar2 = m_betaBar;
  double t1 = (1. + m_betaBar1 * m_betaBar2) / (m_betaBar1 + m_betaBar2);
  double logarg =  2.*(1. + m_betaBar1) * (1. + m_betaBar2) / ((1. - m_betaBar1) * (1. - m_betaBar2));
  m_gBar  = -m_QF2 * m_alpi * t1 * (log(logarg / 2.) - 2.); // See Mareks phd thesis A.2.1
  m_gpBar = -m_QF2 * m_alpi * t1 * (log(logarg / 2.));
  m_fsrform = exp(m_gBar / 4 + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5));
  m_fsrWeight = m_massW * m_hideW * m_wt2 * exp(m_gBar / 4 + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5));
  if (IsBad(m_fsrWeight)) {
    msg_Error() << METHOD << "\n FSR weight is NaN\n"
                << "\n Eprime = " << sqrt(m_dip_sp)
                << "\n Eq = " << sqrt(m_sQ)
                << "\n EminQ = " << m_EminQ
                << "\n q1q2 = " << m_q1q2
                << "\n Exp(YFS) = " << m_expf
                << "\n YFS_IR = " << m_YFS_IR
                << "\n VolMc = " << m_volmc
                << "\n btil = " << m_btil
                << "\n btildestar = " << m_btilStar
                << "\n Mass Weight = " << m_massW
                << "\n dipole = " << m_dipole
                << "\n r1 = " << m_r1
                << "\n r2 = " << m_r2
                << "\n mass r1 = " << m_r1.Mass()
                << "\n mass r2 = " << m_r2.Mass()
                << "\n Hidden Photon Weight = " << m_hideW
                << "\n Photon Scale Weight =  " << m_wt2 << "\n";
  }
  DEBUG_FUNC("FSR for Dipole  = " << m_dipoleFl
             << "\n N Photons = " << m_N
             << "\n N Photons removed = " << m_NRemoved
             << "\n Eprime = " << sqrt(m_dip_sp)
             << "\n Eq = " << sqrt(m_sQ)
             << "\n EminQ = " << m_EminQ
             << "\n q1q2 = " << m_q1q2
             << "\n Exp(YFS) = " << m_expf
             << "\n YFS_IR = " << m_YFS_IR
             << "\n VolMc = " << m_volmc
             << "\n btil = " << m_btil
             << "\n btildestar = " << m_btilStar
             << "\n Mass Weight = " << m_massW
             << "\n dipole = " << m_dipole
             << "\n m_1 = " << m_dipole[0].Mass()
             << "\n m_2 = " << m_dipole[1].Mass()
             << "\n m_v = " << (m_dipole[0] + m_dipole[1]).Mass()
             << "\n r1 = " << m_r1
             << "\n r2 = " << m_r2
             << "\n mass r1 = " << m_r1.Mass()
             << "\n mass r2 = " << m_r2.Mass()
             << "\n Hidden Photon Weight = " << m_hideW
             << "\n Photon Scale Weight =  " << m_wt2
             << "\n Cut is =  " << m_cut
             << "\n Total Weight = " << m_fsrWeight << "\n");
}



void FSR::BoostToXFM() {
  // p_rot   = new Poincare(m_dipole[0],Vec4D(0.,0.,0.,1.));
  Vec4D Q = m_dipole[0] + m_dipole[1] + m_photonSum;
  ATOOLS::Poincare poin(Q);
  ATOOLS::Poincare prot1(m_dipole[0], Vec4D::ZVEC);
  ATOOLS::Poincare prot2(m_dipole[1], -Vec4D::ZVEC);
  Poincare rotate;
  rotate = Poincare(m_dipole[0], Vec4D(0., 0., 0., 1.));
  for (auto &p : m_dipole) {
    poin.Boost(p);
  }
  for (auto &p : m_photons) {
    poin.Boost(p);
  }
  poin.Boost(m_photonSum);
  Q = {sqrt(m_s), 0, 0, 0};
  ATOOLS::Poincare poin2(Q);
  for (auto &p : m_dipole) {
    poin.BoostBack(p);
  }
}

void FSR::RotateDipole() {
  double costh = 1. - 2.*ran->Get();
  double theta = acos(costh);
  double phi = 2.*M_PI * ran->Get();
  Vec4D t1 = m_dipole[0];
  Vec4D t2 = m_dipole[1];
  int i(0);
  Vec4D t;
  for (auto &p : m_dipole) {
    if (i == 0) t = t1;
    else t = t2;
    p[2] = cos(theta) * p[2] - sin(theta) * p[3];
    p[3] = sin(theta) * p[2] + cos(theta) * p[3];

    p[1] = cos(phi) * p[1] - sin(phi) * p[2];
    p[2] = sin(phi) * p[1] + cos(phi) * p[2];
    i++;
  }
  for (auto &p : m_photons) {
    p[2] = cos(theta) * p[2] - sin(theta) * p[3];
    p[3] = sin(theta) * p[2] + cos(theta) * p[3];

    p[1] = cos(phi) * p[1] - sin(phi) * p[2];
    p[2] = sin(phi) * p[1] + cos(phi) * p[2];
  }
  m_photonSum[2] = cos(theta) * m_photonSum[3] - sin(theta) * m_photonSum[3];
  m_photonSum[3] = sin(theta) * m_photonSum[3] + cos(theta) * m_photonSum[3];

  m_photonSum[1] = cos(phi) * m_photonSum[1] - sin(phi) * m_photonSum[2];
  m_photonSum[2] = sin(phi) * m_photonSum[1] + cos(phi) * m_photonSum[2];
}

void FSR::RejectEvent() {
  DEBUG_FUNC("EVENT REJECETED" << " Exp(YFS) = " << m_expf
             << "\n YFS_IR = " << m_YFS_IR
             << "\n VolMc = " << m_volmc
             << "\n Mass Weight = " << m_massW << "\n"
             << "Hidden Photon Weight = " << m_hideW
             << "\n Photon Scale Weight =  " << m_wt2);
  m_f = m_fbar = 0.0;
  m_fsrWeight  = 0.0;
  m_hideW = 0.0;
  m_photonSum *= 0;
  m_photons.clear();
  m_MassWls.clear();
  m_massW = 0.0;
  m_cut = 1;
  m_sprim = 0;
  m_failed = true;
}

void FSR::Reset() {
  m_f = m_fbar = 0.0;
  m_fsrWeight  = 0.0;
  m_hideW = 0.0;
  m_photonSum *= 0;
  m_photons.clear();
  m_MassWls.clear();
  m_massW = 0.0;
  m_cut = 0;
  m_sprim = 0;
  m_cos.clear();
  m_sin.clear();
  m_yini.clear();
  m_zini.clear();
}




double FSR::Eikonal(Vec4D k) {
  return -m_alpi / (4.*M_PI) * (m_dipole[0] / (m_dipole[0] * k) - m_dipole[1] / (m_dipole[1] * k)).Abs2();
}

double FSR::EikonalInterferance(Vec4D k) {
  return m_alpi / (4.*M_PI) * 2.*m_dipole[0] * m_dipole[1] / ((k * m_dipole[0]) * (k * m_dipole[1]));
}
