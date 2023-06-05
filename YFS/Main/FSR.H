#ifndef YFS_Main_FSR_H
#define YFS_Main_FSR_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include "YFS/Main/YFS_Form_Factor.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Tools/Dipole.H"
#include "YFS/Main/Define_Dipoles.H"
// #define FSR_DEBUG

using namespace ATOOLS;

namespace YFS {
  class FSR : public YFS_Base {
    public:
      FSR();
      ~FSR();
      ATOOLS::Poincare *p_rot;
      double m_u, m_fsrcut, m_Kmin, m_Kmax, m_sQ,m_sprim,m_sqtest, m_sX, m_Emin, m_QF2, m_Q1, m_Q2;
      double m_phi, m_theta, m_c, m_eta1, m_eta2, m_EminQ, m_mass1, m_mass2,m_q1q2, m_st, m_p1p2;
      double m_beta, m_beta1, m_beta2, m_mu1, m_mu2, m_g, m_gp, m_EQ, m_betaf, m_omegaF, m_volmc, m_gpBar, m_gBar;
      double m_massW, m_hideW, m_expf, m_wt2, m_amc2, m_amfin, m_EminPP, m_f, m_fbar,m_bvrA,m_pT,m_yy,m_xfact;
    double  m_fsrform, m_nbar, m_betaBar,m_YFS_IR, m_btil, m_btilStar, m_volumc, m_delvol, m_DelYFS,m_delta1,m_cut;
    double m_BtiXcru, m_BtiQcru, m_Edelta, m_fsrcutF, m_omegaMax, m_omegaMin, m_A4, m_A, m_betaBar1, m_betaBar2;
      // double m_btilde
      std::vector<double> m_mass, m_MassWls, m_cos, m_sin, m_dist1, m_dist2, m_del1, m_del2, m_yini, m_zini, m_phi_vec, m_k0, m_rand, m_fvev, m_fbarvec;
    int m_n, m_N, m_NRemoved, m_eikonal_mode, m_use_crude;
      bool m_kkmcAngles, m_use_massive_nbar, m_failed, m_use_model_alpha;
    Vec4D_Vector m_photons, m_dipole, m_photonspreboost;
      Flavour_Vector m_dipoleFl;
      Vec4D m_photonSum, m_photonSumPreBoost, m_r1, m_r2, m_QFrame, m_bornQ1, m_bornQ2, m_Q, m_px;
      YFS::Dipole *p_dipole;
      void GeneratePhotonMomentum();
      void GenerateAngles();
      void DefineGammas();
      double CalculateBeta(const Vec4D& p);
      void CalculateBetaBar();
      bool MakeFSR();
      void Reset();
      void NPhotons();
      double Btild(Vec4D p1,Vec4D p2, double kmax);
      double chi(double eta, double y1, double y2, double y3, double y4);
      double Zij(double eta, double y1, double y2);
      double LI(double x);
      double YEta(double eta,double y1, double y2, double y3, double y4);
      double BVR_A(Vec4D p1, Vec4D p2);
      double Eikonal(Vec4D);
      double EikonalInterferance(Vec4D);
      inline double GetWeight() {return m_fsrWeight;}
      inline double GetSQ() {return m_sQ;}
      inline Vec4D GetPhotonSum() {return m_photonSum; }
      inline Vec4D_Vector GetPhotons() {return m_photons; }
      inline void SetNPhoton(int N) {m_n = N;}
      inline void SetPhotons(Vec4D k) {m_photonSum = k;}
      inline void SetDipole(Vec4D_Vector dip) {m_dipole = dip;}
      inline void SetV(double v) {m_v = v;}
      inline Vec4D_Vector GetDipole() {return m_dipole;}
    YFS::YFS_Form_Factor *p_isrFormFact, *p_fsrFormFact, *p_yfsFormFact;
      void HidePhotons();
      // void MakePair(double cms, Vec4D &p1, Vec4D &p2, double mass1, double mass2);
      void MakePair(double cms, Vec4D &p1, Vec4D &p2);
      void MakePair(double cms, Vec4D &p1, Vec4D &p2, double mass1, double mass2,
                    double &eta1, double &eta2);
      void BoostDipole(Vec4D_Vector &dip);
      void Weight();
      bool F(ATOOLS::Vec4D_Vector &k);
      bool Initialize(Dipole_Vector::iterator dipole);
      void BoostToXFM();
      void RejectEvent();
      void RotateDipole();
      void YFS_FORM();
      void RescalePhotons();
      // PHOTONS::YFS_Form_Factor *p_yfsformfactor;
  };

}

#endif
