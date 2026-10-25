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
      //! Flavours m_flav_label was last built from. m_flav_label is invariant
      //! per process, so MakeDipolesII/MakeDipoles skip rebuilding it (a
      //! std::map, so tree-node allocations) when the flavours have not
      //! changed since last event -- same guard YFS_Handler::SetFlavours uses.
      ATOOLS::Flavour_Vector m_flav_label_fl;
      // std::map<ATOOLS::Vec4D, int> m_mom_label;
      ATOOLS::Vec4D_Vector m_softphotons, m_bornmomenta;
      YFS::DipoleSet m_set;
      std::map<std::string,MODEL::Vertex_List> m_proc_restab_map;
      std::unique_ptr<YFS::YFS_Form_Factor> p_yfsFormFact;
      DivArrD m_virtSub;

      // m_type removed: shadowed YFS_Base's std::string m_type with a
      // different type entirely, and was never read.
      Define_Dipoles();
      ~Define_Dipoles();

      double m_N_born_Gamma;

    
      void AddRadiation();
      void MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born);
      void MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,  ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesFF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,ATOOLS::Vec4D_Vector const &born);
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
      //! the dim.-reg. virtual B of an initial-final dipole (YFS: IF_VIRTUAL_TCHANNEL)
      DivArrD IFVirtualBEps(YFS::Dipole &D);
      //! The real-soft B-tilde of the dipoles CalculateVirtualSub sums, photon
      //! energies below K in the frame of the legs, in the IR regularisation
      //! of that subtraction (dim-reg finite part: lambda^2 = 4 pi mu^2/eps).
      double RealSoftSum(double K);
      //! diagnostic for RealSoftSum: per dipole type and IR-finite flag
      std::string RealSoftSumReport(double K, const Vec4D_Vector &dirs);
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
      //! FormFactor built from an already-computed FormFactorSum(), so a caller
      //! that needs both gets them from one evaluation instead of two.
      double FormFactor(double form);
      double FormFactorSum();
      //! FormFactorSum() without its initial-state (II) term: the final-
      //! state, Coulomb, initial-final and decay pieces, on the legs the
      //! dipoles were built with (RV_MODE remainder_own_legs)
      double FormFactorSumFinalLegs();
      //! Just the initial-final term of FormFactorSum(), i.e. the piece
      //! IFI_Sub=1 adds and IFI_Sub=0 leaves out. Zero when IFI_Sub is off.
      //! Lets the IFI rung be undone exactly as exp(-FormFactorSumIF())
      //! instead of costing a separate production run.
      double FormFactorSumIF();
      //! The same sum with the soft cutoff at omega instead of IFIOmega().
      double FormFactorSumIF(double omega);
      //! The initial and final real B-tilde terms of FormFactorSum() moved
      //! from their exponent regions to the common generation region (photon
      //! energy below emin in the beams' rest frame). ACRAIC: FIXED_ORDER_STREAM
      double FormFactorShiftToRegion(double emin);
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
      //! The W's and their daughters' positions, as BuildPole found them.
      inline const DipoleSet::WWLegs &WW() const { return m_set.WW(); }
      /*!
        Once per run: the closed-form decay exponent (YFS_Form_Factor::
        BVR_decay, what FormFactorSumDecay exponentiates) against the same
        dipole's virtual plus real B from the dim-reg machinery in the W rest
        frame - the pieces the O(alpha) subtraction uses. They must agree, or
        the expansion of exp(Y) and the subtracted V and R are not the same
        object at O(alpha).
      */
      void CheckDecayFormFactor(YFS::Dipole &D);
      bool m_decaycheckdone = false;
      //! YFS: WW_DECAY_VIRTUAL_SUB, see the definition
      bool DecayVirtualSubtraction() const;
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
      inline bool HasDipoleII() const { return m_set.HasII(); }
      inline Dipole & GetDipoleII() { return m_set.II(); }
      inline YFS::DipoleView GetDipoleFF() { return m_set.FF(); }
      inline void SetProcResMap(std::map<std::string,MODEL::Vertex_List> map) {m_proc_restab_map=map;}
      inline double Get_Finite() {return m_virtSub.GetFinite();}
      inline double Get_E1() {return m_virtSub.GetIR();}
  };


}
#endif
