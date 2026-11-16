#ifndef YFS_Main_Define_Dipoles_H
#define YFS_Main_Define_Dipoles_H

#include "ATOOLS/Phys/Particle.H"
#include "YFS/Main/Dipole.H"
#include <memory>
#include "YFS/Main/DipoleSet.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/YFS_Form_Factor.H"


namespace MODEL {
  class Single_Vertex;
  typedef std::vector<Single_Vertex *> Vertex_List;
}

namespace YFS {
  typedef std::vector<YFS::Dipole> Dipole_Vector;


  class Define_Dipoles : YFS_Base {
    public:
    void SetS(const double s) override;
      int m_N;
      std::map<ATOOLS::Flavour, int> m_flav_label;
      //! Flavours m_flav_label was last built from (rebuilt only on change).
      ATOOLS::Flavour_Vector m_flav_label_fl;
      // std::map<ATOOLS::Vec4D, int> m_mom_label;
      ATOOLS::Vec4D_Vector m_softphotons, m_bornmomenta;
      YFS::DipoleSet m_set;
      std::map<std::string,MODEL::Vertex_List> m_proc_restab_map;
      std::unique_ptr<YFS::YFS_Form_Factor> p_yfsFormFact;
      DivArrD m_virtSub;

      Define_Dipoles();
      ~Define_Dipoles();

      double m_N_born_Gamma;

    
      void AddRadiation();
      void MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born);
      void MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,  ATOOLS::Vec4D_Vector const &born);
      void Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom);
      // Flags the subset of m_dipolesFF that may radiate: a matching of the
      // charged final state, so each leg enters FSR exactly once.
      double ResonanceWidthDistance(YFS::Dipole &D);
      double ResonantDist(YFS::Dipole &D, const Vec4D &k);
      double ResonantDist(YFS::Dipole &D);
      dipoletype::code WhichResonant(const Vec4D &k);
      double CalculateRealSub(const Vec4D &k);
      double CalculateRealSubIF(const Vec4D &k);
      double CalculateRealSubEEX(const Vec4D &k);
      double CalculateVirtualSub();
      //! diagnostic (YFS: ACRAIC_Compare): CalculateVirtualSub's default branch
      //! split by dipole type, initial / final / initial-final
      void VirtualSubParts(double &ii, double &ff, double &ifi);
      //! the dim.-reg. virtual B of an initial-final dipole (YFS: IF_VIRTUAL)
      DivArrD IFVirtualBEps(YFS::Dipole &D);
      //! The real-soft B-tilde of the dipoles CalculateVirtualSub sums, photon
      //! energies below K in the frame of the legs, in the IR regularisation
      //! of that subtraction (dim-reg finite part: lambda^2 = 4 pi mu^2/eps).
      double RealSoftSum(double K);
      double CalculateVirtualSubEps();
      double CalculateRealVirtualSubEps(const Vec4D &k);
      double CalculateEEX();
      double CalculateEEXVirtual();
      double CalculateEEXVirtual(int betaorder);
      double CalculateEEX(const Vec4D & k1, const Vec4D & k2);
      double CalculateFlux(const Vec4D &k);
      double CalculateFlux(const Vec4D &k, dipoletype::code &fluxtype);
      double CalculateFlux(const Vec4D &k, const Vec4D &kk);
      double CalculateFlux(const Vec4D &k, const int &fsrmode);
      double Propagator(const double &s, int width=0);
      double FormFactor();
      //! FormFactor from an already-computed FormFactorSum()
      double FormFactor(double form);
      double FormFactorSum();
      //! FormFactorSum() without its initial-state (II) term, on the legs the
      //! dipoles were built with (RV_MODE remainder_own_legs)
      double FormFactorSumFinalLegs();
      //! The initial-final term of FormFactorSum() alone (zero when IFI_Sub is off)
      double FormFactorSumIF();
      //! The same sum with the soft cutoff at omega instead of IFIOmega().
      double FormFactorSumIF(double omega);
      //! The initial and final real B-tilde terms of FormFactorSum() moved
      //! from their exponent regions to the common generation region (photon
      //! energy below emin in the beams' rest frame). ACRAIC: FIXED_ORDER_STREAM
      double FormFactorShiftToRegion(double emin);
      //! the final-state dipoles' part of FormFactorShiftToRegion (ACRAIC: FIXED_ORDER_EVENT_EXPONENT)
      double FormFactorShiftToRegionFinal(double emin);
      //! The decay-stage term of the WW pole expansion (zero in the flat scheme)
      double FormFactorSumDecay();
      //! Coulomb singularity removed from a massive-vector-pair dipole's form
      //! factor, eq (9) of hep-ph/9606429. Zero unless the pole expansion is
      //! active and the pair is slow.
      double CoulombSubtraction(YFS::Dipole &D);
      //! Pole-expansion attempts and events inside the leading-pole window
      long m_pole_tried = 0, m_pole_built = 0;
      bool m_pole_active = false;
      long m_recoil_tried = 0, m_recoil_massfail = 0;
      //! Whether the pole expansion is live for the current event.
      inline bool PoleActive() const { return m_pole_active; }
      //! The W's and their daughters' positions, as BuildPole found them.
      inline const DipoleSet::WWLegs &WW() const { return m_set.WW(); }
      /*!
        Once per run: compare the closed-form decay exponent (BVR_decay)
        with the dipole's dim-reg virtual plus real B in the W rest frame.
      */
      void CheckDecayFormFactor(YFS::Dipole &D);
      bool m_decaycheckdone = false;
      //! YFS: WW_DECAY_VIRTUAL_SUB, see the definition
      bool DecayVirtualSubtraction() const;
      //! Carry the production-stage recoil from the W's down to the four
      //! fermions. Call after the W-pair dipole has radiated.
      bool ApplyPoleRecoil(ATOOLS::Vec4D_Vector &plab);
      bool MakeDipolesPole(ATOOLS::Flavour_Vector const &fl,
                           ATOOLS::Vec4D_Vector const &mom,
                           ATOOLS::Vec4D_Vector const &born);
      //! One-shot dump (first call only) of the constructed dipole set:
      //! flavours, charges, radiating pairs, soft cutoffs. YFS: Dump_Dipoles.
      void DumpDipoles();
      double IFIOmega() const;
      double RealIFWeight(const ATOOLS::Vec4D_Vector &photons);
      // YFS_Base is inherited privately, so expose what the handler reports.
      inline long   IFIClipped() const { return m_ifi_clipped; }
      inline double IFIRClip()   const { return m_ifi_rclip; }
      inline bool HasDipoleII() const { return m_set.HasII(); }
      inline Dipole & GetDipoleII() { return m_set.II(); }
      inline YFS::DipoleView GetDipoleFF() { return m_set.FF(); }
      inline void SetProcResMap(std::map<std::string,MODEL::Vertex_List> map) {m_proc_restab_map=map;}
      inline double Get_E1() {return m_virtSub.GetIR();}
  };


}
#endif
