#ifndef METOOLS_Current_Library_VA_P_A_H
#define METOOLS_Current_Library_VA_P_A_H

#include "METOOLS/HadronCurrents/VA_P_X_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// D -> A l nu, X with J^P = 1^+  (K_1(1270), K_1(1400)).
  /// charm_semileptonic.tex Sec."D to A", Eqs.(PAvector),(PAaxial):
  /// the parity-flipped analogue of VA_P_V, with {V_0,V_1,V_2} in the
  /// vector current and a single A in the axial one, and with the
  /// roles of (M+m) and (M-m) interchanged relative to D->V.
  ///
  /// As in VA_P_V the q^mu pieces are collected before implementation,
  ///
  ///   V^mu = -i { (eps*.q) G(q^2)/q^2 q^mu
  ///               + (M-m) V_1(q^2) eps*^mu
  ///               - V_2(q^2)(eps*.q)/(M-m) P^mu },
  ///   G(q^2) = 2m V_0 - (M-m)V_1 + (M+m)V_2 = 2m[V_0-V_3],
  ///
  /// using Delta/(M-m)=M+m, so that Eq.(V0constraint) V_0(0)=V_3(0)
  /// makes G(0)=0 exactly. V_0(0) is computed per event, never read.
  ///
  /// Normalisation convention follows the BESIII K_1(1270) dynamics
  /// fit, which measures only the RATIOS
  ///     rA   = A(0)/V_1(0)   = -0.112 +- 0.010 +- 0.009,
  ///     rV2  = V_2(0)/V_1(0) = -0.043 +- 0.010 +- 0.025,
  /// with no absolute V_1(0). Since the channel rate is fixed by the
  /// branching ratio in the decay table anyway, only the ratios affect
  /// anything observable here; V1_0 defaults to 1.
  ///
  /// !! SIGN WARNING !! the note flags explicitly that the signs of
  /// those two ratios may only be transported after matching the
  /// precise D->A vector/axial and state-phase conventions. A global
  /// rephasing of the K_1 is harmless; a relative vector/axial sign is
  /// not. Use SIGN_A / EPS_SIGN to flip if the closure test demands it.
  //////////////////////////////////////////////////////////////////////////
  class VA_P_A : public VA_P_X_Base {
  private:
    FF_P_X * p_A, * p_V0, * p_V1, * p_V2;
    double   m_V1_0, m_rA, m_rV2;
    bool     m_useratios;
    /// Overall phase of the vector piece, exp(i*m_phase); default
    /// -pi/2, i.e. the "-i" of Eq.(PAvector).
    double   m_phase;
  public:
    VA_P_A(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name);
    ~VA_P_A() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
