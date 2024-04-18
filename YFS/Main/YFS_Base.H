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

  std::istream &operator>>(std::istream &str,submode::code &sm);

  class YFS_Base { 
  protected:

    double m_exponent, m_xmax;
    bool   m_on;

  public:
    YFS_Base();
    virtual ~YFS_Base();

    bool On() const   { return m_on;   }

    double m_mass2[2];
    double m_smin;
    double m_s;
    double m_v, m_vmin, m_vmax, m_vmin1, m_isrcut, m_sp, m_min_photon_E;
    double m_alpha, m_alpi, m_hardmin,m_gp,m_g, m_pole_fac;
    bool m_CalForm, m_bornInit, m_fillblob, m_isr_debug, m_fsr_debug, m_born_set, m_coulomb, m_constfsrW, m_setparticles;
    double m_photonMass, m_oneloop, m_born,m_isrWeight, m_fsrWeight, m_rescale_alpha;
    double m_nbar, m_deltacut, m_betatWW, m_eps, m_fixed_weight, m_resonace_max;
    int m_yfsmode, m_N, m_order, m_nmax, m_betaorder, m_N_FSR, m_fullform, m_formWW, m_ifisub;
    int m_fsrmode, m_useint, m_check_mass_reg, m_check_poles, m_check_real,m_check_rr, m_use_fsr_beta;
    int m_mode, m_rmode, m_semiyfs, m_NReal, m_qedmode,m_kkmcAngles, m_griff, m_fixed_photons; 
    int m_no_born, m_check_real_sub, m_no_subtraction,m_tchannel, m_resonance_mode, m_check_virt_born,m_noflux, m_massless_sub;
    bool m_looptool, m_setmass, m_coll_real, m_RealPhotons, m_virtual_only, m_real_only, m_use_model_alpha;
    bool m_beamspread, m_isrinital,m_realtool,m_realvirt,m_useceex,m_realrealtool, m_check_rr_sub;
    int m_flux_mode, m_hidephotons, m_photon_split, m_int_nlo;
    std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR, m_debugDIR_NLO;
    Vec4D m_beam1, m_beam2;
    // beam energy spread values
    double m_beam_sigx, m_beam_sigy, m_beam_rho;
    submode::code m_submode;


    void RegisterSettings();
    void RegisterDefaults();
    double Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
    double EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
    // YFS::Coulomb *p_coulomb;
    // PHASIC::Process_Base *proc;
    // ATOOLS::Particle_Vector m_particles,  dipoleISR, dipoleFSR;
    // ATOOLS::ParticleMomMap m_inparticles, m_outparticles;

    // std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR;
    // ATOOLS::Poincare m_cmsboost;
    // ATOOLS::Vec4D    m_fixvecs[2], p_cms[2], m_realPho, m_Q,m_P, m_PX, m_Wm, m_Wp;
    // ATOOLS::Flavour m_bunch;
    // ATOOLS::Flavour_Vector m_flavs;
    // std::vector<double> m_info_lab, m_info_cms, m_dvir, m_y, m_z, m_eikonalsISR, m_eikonalsFSR;
    // Remnant_Base    *p_remnants[2];

    double m_f;

  };
}

#endif
