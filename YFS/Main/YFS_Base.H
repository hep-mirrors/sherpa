#ifndef YFS_Main_YFS_Base_H
#define YFS_Main_YFS_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"


using namespace ATOOLS;

namespace YFS {
struct submode {
  enum code {
    off        = 0,
    local       = 1,
    global       = 2
  };
};
struct yfsmode {
  enum code {
    off       = 0,
    isr       = 1,
    isrfsr    = 2,
    fsr       = 3,
  };
};

struct wgt {
  enum code {
    off       = 0,
    full      = 1,
    mass      = 2,
    hide      = 3,
    jacob     = 4,
  };
};

struct fixed_order {
  enum code {
    off       = 0,
    full      = 1,
    lo        = 2,
    nlo       = 3,
    nnlo      = 4,
  };
};



std::istream &operator>>(std::istream &str, submode::code &sm);
std::istream &operator>>(std::istream &str, yfsmode::code &sm);
std::istream &operator>>(std::istream &str, wgt::code &sm);
std::istream &operator>>(std::istream &str, fixed_order::code &sm);
std::ostream &operator<<(std::ostream &str, const yfsmode::code &sm);
std::ostream &operator<<(std::ostream &str, const wgt::code &sm);
std::ostream &operator<<(std::ostream &str, const fixed_order::code &sm);

class YFS_Base {
protected:

  double m_exponent, m_xmax;
  bool   m_on;

public:
  YFS_Base();
  virtual ~YFS_Base();

  bool On() const   { return m_on;   }

  double m_mass2[2];
  double m_s,m_irscale;
  double m_v, m_vmin, m_vmax, m_isrcut, m_sp;
  double m_alpha, m_alpi, m_hardmin, m_gp, m_g;
  bool m_CalForm, m_fillblob, m_isr_debug, m_fsr_debug, m_coulomb, m_setparticles;
  double m_photonMass, m_oneloop, m_born, m_isrWeight, m_fsrWeight, m_rescale_alpha;
  double m_nbar, m_deltacut, m_betatWW, m_eps, m_fixed_weight, m_resonace_max;
  double m_rv_test_x, m_rv_test_theta, m_rv_test_phi;
  double m_rv_soft_cut; // soft-photon energy fraction below which RV is skipped
  // Relative-cancellation threshold for the RV beta_1^1. The RV matrix element
  // m_rv and its subtraction m_rvsub each diverge like S(k)~1/k^2 and cancel to
  // form d1; when |d1| drops below m_rv_cancel_eps*max(|m_rv|,|m_rvsub|) the
  // result is roundoff-dominated (and, being that small, physically
  // negligible), so the photon is skipped. <=0 disables. This is the
  // process-independent replacement for the flat RV_SOFT_CUT energy gate.
  double m_rv_cancel_eps;
  // Matrix-element stability guard for the RV. The one-loop provider is
  // numerically unstable for soft, beam-collinear (ISR) photons, returning a
  // real-virtual ME rv up to ~1e30x larger than the scale its own subtraction
  // predicts (max(|rvsub|,|aB|)); the failures are confined to soft/collinear
  // (verified: hard wide-angle photons are always stable), so skipping them
  // loses negligible physics. Skip when |rv| > m_rv_me_max_ratio*subscale.
  // Recommended O(1e6) - the healthy/unstable ratio gap spans 1e4..1e25.
  // <=0 disables.
  double m_rv_me_max_ratio;

  int  m_order, m_betaorder, m_fullform, m_formWW, m_ifisub;
  // When set, fill the in-run RV plateau-scan histograms (RV contribution
  // binned by log10 cancellation ratio and by photon energy fraction) so a
  // single production run yields RV_total(cut) and the cut can be validated.
  int  m_rv_cancel_hist;
  int  m_check_mass_reg, m_check_poles, m_check_real, m_check_rr, m_check_rv, m_rv_hard_photon;
  int m_rmode, m_qedmode, m_kkmcAngles, m_virt_sub,m_mingammaN;
  int m_no_born, m_check_real_sub, m_check_rr_sub, m_no_subtraction, m_tchannel, m_check_virt_born, m_noflux, m_massless_sub;
  
  bool m_looptool, m_virtual_only, m_real_only, m_use_model_alpha,m_failcut,m_nlocuts;
  bool m_isrinital, m_realtool, m_realvirt, m_useceex, m_coll_real, m_rrtool, m_vvtool, m_dim_reg;
  bool m_skipNegWeights;
  bool m_fb_analysis;
  int m_fb_kf;
  int m_flux_mode, m_hidephotons, m_photon_split, m_int_nlo, m_eex_virt;
  std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR, m_debugDIR_NLO;
  Vec4D m_beam1, m_beam2;

  submode::code m_submode;
  yfsmode::code m_mode;
  fixed_order::code m_fixedOrder;


  void RegisterSettings();
  void RegisterDefaults();
  double Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
  double EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
  inline bool  HasISR()   const {return (m_mode != yfsmode::off && m_mode != yfsmode::fsr); }
  inline bool  HasFSR()   const {return (m_mode != yfsmode::off && m_mode != yfsmode::isr); }
  inline fixed_order::code FixedOrder()  const {return m_fixedOrder;}

};
}

#endif
