#ifndef YFS_Main_Generate_Photon_H
#define YFS_Main_Generate_Photon_H

#include "ATOOLS/Phys/Particle.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Histogram.H"

#include "YFS/Main/YFS_Base.H"
#include "YFS/Tools/Dipole.H"

using namespace ATOOLS;
namespace BEAM   { class Beam_Base;      }
// #define YFS_DEBUG

namespace YFS {

#ifdef YFS_DEBUG
typedef std::map<std::string, ATOOLS::Histogram *> histograms;
static histograms m_histograms;
#endif
class ISR : public YFS_Base {
public:
  double  m_b1, m_sp;
  double  m_b2;
  double  m_theta, m_angleWeight;
  double  m_phi;
  double  m_c, m_F, m_sin;
  double  m_w, m_f, m_fbar;
  double  m_delta, m_betaf;
  double  m_jacW, m_Jac, m_Jac0, m_lam, m_lam0, m_diljac, m_diljac0;
  double  A;
  int m_n, m_cut;
  double m_beta, m_betap, m_form, m_formF, m_weight, m_Kmax, m_Kmin, m_betak;
  double  m_L, m_am2, m_mass, m_mass2, m_x;
  double m_massW;
  double m_W, m_W0; // for beta1
  bool m_setmass;
  double m_nsuccess, m_nfail, m_ntotal;
  // ATOOLS::Particle *p_photon;
  std::vector<double> m_yini, m_AA, m_K2, m_PTK, m_del1, m_del2;
  std::vector<double> m_zini, m_eikonals, m_cos, m_jacvec, m_scale;
  Vec4D m_pvec1, m_pvec2, m_qvec1, m_qvec2, m_photon, m_photonSum;
  Vec4D_Vector m_photons;
  void GeneratePhotonAngleMassless();
  BEAM::Beam_Base **p_beam;


public:
  ISR();
  ~ISR();
  void GeneratePhotonMomentum();
  double Jacobian( double x);
  double GenerateV();
  void SetIncoming(YFS::Dipole *p_dipole);
  double CalculateBeta(const ATOOLS::Vec4D &p);
  void GenerateAngles();
  void Weight();
  void WeightFSR();
  double EikonalWeight(const Vec4D_Vector, const Vec4D_Vector, double a);
  double Eikonal(Vec4D);
  void NPhotons();
  void MakeYFS();
  void Beta1();
  void D2(Vec4D, double, double, double);
  void MapPhotonMomentun();
  void MapMomenta(Vec4D_Vector p, Vec4D_Vector q, Vec4D_Vector k);
  void Sort(Vec4D_Vector &p);
  void Sort(std::vector<double> &p);
  void GenerateAnglesForMassless();
  double Flux();
  double Eikonal(Vec4D k, Vec4D p1, Vec4D p2);
  void Clean();

  inline double GetWeight() { return m_weight; }
  inline ATOOLS::Vec4D GetPhotonSum() { return m_photonSum; }
  inline ATOOLS::Vec4D_Vector GetPhotons() { return m_photons; }
  inline double  GetTheta() { return m_theta; }
  inline double  GetPhi() { return m_phi; }
  inline double  GetN() { return m_n; }
  inline void  SetV(double v) { m_v = v; }
  inline void  Setmass(bool setmass) { m_setmass = setmass; }
  inline double  GetGamma() { return m_g; }
  inline double  GetGammaPrime() { return m_gp; }
  inline int  GetCut() { return m_cut; }
  inline std::vector<double>  GetSudakovY() { return m_yini; }
  inline std::vector<double>  GetSudakovZ() { return m_zini; }
  inline std::vector<double>  GetEikonals() { return m_eikonals; }
  void SetMinimalPhotonAngle(const double&);


};



}

#endif
