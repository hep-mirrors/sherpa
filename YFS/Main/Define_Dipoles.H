#ifndef YFS_Main_Define_Dipoles_H
#define YFS_Main_Define_Dipoles_H

#include "ATOOLS/Phys/Particle.H"
#include "YFS/Main/Dipole.H"
#include "YFS/Main/DipoleSet.H"
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
    
    // Individual dipole creation
    void CreateFinalStateDipole(const ParticleInfo& particle1,
                               const ParticleInfo& particle2);
    void CreateResonantDipole(const ParticleInfo& particle1,
                             const ParticleInfo& particle2);
    
    // State management
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
      YFS::DipoleSet m_set;
      std::vector<std::vector<int>> m_dip;
      std::map<std::string,MODEL::Vertex_List> m_proc_restab_map;
      YFS::YFS_Form_Factor *p_yfsFormFact;
      DivArrD m_virtSub, m_vvSub;

      // m_type removed: shadowed YFS_Base's std::string m_type with a
      // different type entirely, and was never read.
      Define_Dipoles();
      ~Define_Dipoles();

      int m_in, m_out;
      double m_N_born_Gamma;
      std::vector<double> m_massInC, m_massInN;
      std::vector<double> m_massOutC, m_massOutN;

      ATOOLS::Flavour_Vector m_dipoleFlav;
      ATOOLS::Vec4D_Vector m_dipoleMom;
    
      void AddRadiation();
      void MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born);
      void MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,  ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesFF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,ATOOLS::Vec4D_Vector const &born);
      void Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom);
      void CleanInParticles();
      void CleanOutParticles();
      // void Get4Mom(ATOOLS::Flavour_Vector const &fl,ATOOLS::Vec4D_Vector mom, std::vector<int> v);
      void Get4Mom(ATOOLS::Flavour_Vector const &fl,ATOOLS::Vec4D_Vector mom);
      bool CheckMasses();
      void IsResonant(YFS::Dipole &D);
      bool IsResonant();
      // Flags the subset of m_dipolesFF that may radiate: a matching of the
      // charged final state, so each leg enters FSR exactly once.
      double ResonanceWidthDistance(YFS::Dipole &D);
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
      double CalculateEEXVirtual(int betaorder);
      double CalculateEEX(const Vec4D & k1, const Vec4D & k2);
      double CalculateFlux(const Vec4D &k);
      double CalculateFlux(const Vec4D &k, dipoletype::code &fluxtype);
      double CalculateFlux(const Vec4D &k, const Vec4D &kk);
      double CalculateFlux(const Vec4D &k, const int &fsrmode);
      double Propagator(const double &s, int width=0);
      double CalculateVirtualSubTchannel();
      double CalculateVirtualSubTchannelEps();
      double FormFactor();
      double FormFactorSum();
      //! Just the initial-final term of FormFactorSum(), i.e. the piece
      //! IFI_Sub=1 adds and IFI_Sub=0 leaves out. Zero when IFI_Sub is off.
      //! Lets the IFI rung be undone exactly as exp(-FormFactorSumIF())
      //! instead of costing a separate production run.
      double FormFactorSumIF();
      //! The decay-stage term of the WW pole expansion. Zero in the flat
      //! scheme.
      double FormFactorSumDecay();
      //! The Coulomb singularity removed from a massive-vector-pair dipole's
      //! form factor, eq (9) of hep-ph/9606429. Zero unless the pole
      //! expansion is active and the pair is slow.
      double CoulombSubtraction(YFS::Dipole &D);
      //! Replace the final-state and interference dipoles with the WW pole
      //! expansion. Returns whether it applied.
      //! How often the pole expansion was attempted and how often the event
      //! was inside the leading-pole window. Reported at the end of the run.
      long m_pole_tried = 0, m_pole_built = 0;
      bool m_pole_active = false;
      long m_recoil_tried = 0, m_recoil_massfail = 0;
      //! Whether the pole expansion is live for the current event.
      inline bool PoleActive() const { return m_pole_active; }
      //! Carry the production-stage recoil from the W's down to the four
      //! fermions. Must be called after the W-pair dipole has radiated and
      //! instead of writing its new momenta straight into the event, since the
      //! W's are not entries in the event record.
      bool ApplyPoleRecoil(ATOOLS::Vec4D_Vector &plab);
      bool MakeDipolesPole(ATOOLS::Flavour_Vector const &fl,
                           ATOOLS::Vec4D_Vector const &mom,
                           ATOOLS::Vec4D_Vector const &born);
      //! One-shot dump of the dipole set as it was actually constructed:
      //! flavours, charges, which pairs radiate, and the soft cutoff each
      //! form factor uses. Prints on the first call and never again. Gated by
      //! the caller on YFS: Dump_Dipoles.
      void DumpDipoles();
      double IFIOmega() const;
      double RealIFWeight(const ATOOLS::Vec4D_Vector &photons);
      // YFS_Base is inherited privately, so expose what the handler reports.
      inline long   IFIClipped() const { return m_ifi_clipped; }
      inline double IFIRClip()   const { return m_ifi_rclip; }
      double TFormFactor();
      double EEXRealVirtual(const Vec4D &k);
      void generate_pairings(std::vector<std::vector<int>>& pairings, std::vector<int>& curr_pairing, std::vector<int>& available_nums);
      void CleanUp();
      std::ostream& operator<<(std::ostream &out);
      inline bool DoneSuccessfully()     { return m_success; }
      inline bool AddedAnything()        { return m_photonsadded; }
      inline void AddSoftPhotons(ATOOLS::Vec4D k) { m_softphotonSum +=k; }
      inline bool HasDipoleII() const { return m_set.HasII(); }
      inline Dipole & GetDipoleII() { return m_set.II(); }
      inline YFS::DipoleView GetDipoleFF() { return m_set.FF(); }
      inline YFS::DipoleView GetDipoleIF() { return m_set.IF(); }
      inline void SetProcResMap(std::map<std::string,MODEL::Vertex_List> map) {m_proc_restab_map=map;}
      inline double Get_Finite() {return m_virtSub.GetFinite();}
      inline double Get_E1() {return m_virtSub.GetIR();}
      inline double GetVV_E1() {return m_vvSub.GetIR();}
      inline double GetVV_E2() {return m_vvSub.GetIR2();}
  };


}
#endif
