#ifndef YFS_Main_Define_Dipoles_H
#define YFS_Main_Define_Dipoles_H

#include "ATOOLS/Phys/Particle.H"
#include "YFS/Tools/Dipole.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/YFS_Form_Factor.H"


namespace MODEL {
  class Single_Vertex;
  typedef std::vector<Single_Vertex *> Vertex_List;
}

namespace {
    constexpr size_t INITIAL_STATE_PARTICLES = 2;
    constexpr size_t TWO_BODY_FINAL_STATE = 4;
}

// Particle information container
struct ParticleInfo {
    ATOOLS::Flavour flavor;
    ATOOLS::Vec4D momentum;
    ATOOLS::Vec4D born_momentum;
    size_t index;
    bool is_initial_state;
    
    ParticleInfo(const ATOOLS::Flavour& f, 
                 const ATOOLS::Vec4D& p, 
                 const ATOOLS::Vec4D& b, 
                 size_t i, 
                 bool initial = false)
        : flavor(f), momentum(p), born_momentum(b), 
          index(i), is_initial_state(initial) {}
};

namespace YFS {
  typedef std::vector<YFS::Dipole> Dipole_Vector;


  class Define_Dipoles : YFS_Base {
    public:
    std::vector<ParticleInfo> ExtractChargedParticles(
        const ATOOLS::Flavour_Vector& flavors,
        const ATOOLS::Vec4D_Vector& momenta,
        const ATOOLS::Vec4D_Vector& born_momenta,
        size_t start_index,
        size_t end_index,
        bool is_initial_state = false) const;
    
    std::vector<ParticleInfo> ExtractInitialStateCharged(
        const ATOOLS::Flavour_Vector& flavors,
        const ATOOLS::Vec4D_Vector& momenta,
        const ATOOLS::Vec4D_Vector& born_momenta) const;
    
    std::vector<ParticleInfo> ExtractFinalStateCharged(
        const ATOOLS::Flavour_Vector& flavors,
        const ATOOLS::Vec4D_Vector& momenta,
        const ATOOLS::Vec4D_Vector& born_momenta) const;
    
    Dipole CreateDipole(
        const ParticleInfo& particle1,
        const ParticleInfo& particle2,
        dipoletype::code type) const;
    
    // Dipole creation by type
    void CreateInitialDipoles(const std::vector<ParticleInfo>& initial_particles);
    void CreateFinalDipoles(const std::vector<ParticleInfo>& final_particles);
    void CreateInitialFinalDipoles(const std::vector<ParticleInfo>& initial_particles,
                                   const std::vector<ParticleInfo>& final_particles);
    
    // Individual dipole creation
    void CreateFinalStateDipole(const ParticleInfo& particle1,
                               const ParticleInfo& particle2);
    void CreateResonantDipole(const ParticleInfo& particle1,
                             const ParticleInfo& particle2);
    
    // State management
    void ResetAllDipoleState();
    // public:
      bool                    m_success;
      bool                    m_photonsadded;
      int m_N;
      std::map<ATOOLS::Flavour, ATOOLS::Vec4D> m_test_dip, m_born_dip;
      std::map<ATOOLS::Flavour, int> m_flav_label;
      // std::map<ATOOLS::Vec4D, int> m_mom_label;
      ATOOLS::Vec4D_Vector m_chargedinparticles;
      ATOOLS::Vec4D_Vector m_neutralinparticles;
      ATOOLS::Vec4D_Vector m_chargedoutparticles;
      ATOOLS::Vec4D_Vector m_neutraloutparticles;
      ATOOLS::Vec4D_Vector m_softphotons, m_bornmomenta;
      ATOOLS::Vec4D m_softphotonSum;
      Dipole_Vector m_dipoles, m_olddipoles, m_newdipoles, m_dipolesII, m_dipolesFF, m_dipolesIF;
      std::vector<std::vector<int>> m_dip;
      std::map<std::string,MODEL::Vertex_List> m_proc_restab_map;
      YFS::YFS_Form_Factor *p_yfsFormFact;
      DivArrD m_virtSub, m_vvSub;

      YFS::dipoletype     m_type;
      Define_Dipoles();
      ~Define_Dipoles();

      int m_in, m_out;
      double m_N_born_Gamma;
      std::vector<double> m_massInC, m_massInN;
      std::vector<double> m_massOutC, m_massOutN;

      ATOOLS::Flavour_Vector m_dipoleFlav;
      ATOOLS::Vec4D_Vector m_dipoleMom;
      void CreateAllDipoles(const ATOOLS::Flavour_Vector& flavors,
                         const ATOOLS::Vec4D_Vector& momenta,
                         const ATOOLS::Vec4D_Vector& born_momenta);
    
      void AddRadiation();
      void MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born);
      void MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,  ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesFF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,ATOOLS::Vec4D_Vector const &born);
      void Dipole_II(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom);
      void Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom);
      void Dipole_IF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void CleanInParticles();
      void CleanOutParticles();
      // void Get4Mom(ATOOLS::Flavour_Vector const &fl,ATOOLS::Vec4D_Vector mom, std::vector<int> v);
      void Get4Mom(ATOOLS::Flavour_Vector const &fl,ATOOLS::Vec4D_Vector mom);
      bool CheckMasses();
      void IsResonant(YFS::Dipole &D);
      bool IsResonant();
      bool CheckResonant(YFS::Dipole &D);
      bool CheckResonant();
      double ResonantDist(YFS::Dipole &D, const Vec4D &k);
      double ResonantDist(YFS::Dipole &D);
      dipoletype::code WhichResonant(const Vec4D &k);
      double CalculateRealSub(const Vec4D &k);
      double CalculateRealSubIF(const Vec4D &k);
      double CalculateRealSubEEX(const Vec4D &k);
      double CalculateVirtualSub();
      double CalculateVirtualSubEps();
      double CalculateVVSubEps();
      double CalculateRealVirtualSubEps(const Vec4D &k);
      double CalculateRealVirtualSub(const Vec4D & k);
      double CalculateEEX();
      double CalculateEEXVirtual();
      double CalculateEEX(const Vec4D & k1, const Vec4D & k2);
      double CalculateFlux(const Vec4D &k);
      double CalculateFlux(const Vec4D &k, dipoletype::code &fluxtype);
      double CalculateFlux(const Vec4D &k, const Vec4D &kk);
      double CalculateFlux(const Vec4D &k, const int &fsrmode);
      double Propagator(const double &s, int width=0);
      double CalculateVirtualSubTchannel();
      double CalculateVirtualSubTchannelEps();
      double FormFactor();
      double TFormFactor();
      double EEXRealVirtual(const Vec4D &k);
      void generate_pairings(std::vector<std::vector<int>>& pairings, std::vector<int>& curr_pairing, std::vector<int>& available_nums);
      void CleanUp();
      std::ostream& operator<<(std::ostream &out);
      inline bool DoneSuccessfully()     { return m_success; }
      inline bool AddedAnything()        { return m_photonsadded; }
      inline void AddSoftPhotons(ATOOLS::Vec4D k) { m_softphotonSum +=k; }
      inline Dipole_Vector Old() { return m_olddipoles; }
      inline Dipole_Vector New() { return m_newdipoles; }
      inline Dipole * GetDipoleII() { return &m_dipolesII[0]; }
      inline Dipole_Vector * GetDipoleFF() { return &m_dipolesFF; }
      inline Dipole_Vector * GetDipoleIF() { return &m_dipolesIF; }
      inline void SetProcResMap(std::map<std::string,MODEL::Vertex_List> map) {m_proc_restab_map=map;}
      inline double Get_Finite() {return m_virtSub.GetFinite();}
      inline double Get_E1() {return m_virtSub.GetIR();}
      inline double GetVV_E1() {return m_vvSub.GetIR();}
      inline double GetVV_E2() {return m_vvSub.GetIR2();}
  };


}
#endif
