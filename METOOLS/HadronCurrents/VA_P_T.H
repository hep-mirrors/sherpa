#ifndef METOOLS_Current_Library_VA_P_T_H
#define METOOLS_Current_Library_VA_P_T_H

#include "METOOLS/HadronCurrents/VA_P_X_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// D -> T l nu, X with J^P = 2^+  (K*_2(1430)).
  /// charm_semileptonic.tex Sec."D to T": exactly the D->V Lorentz
  /// basis after the substitutions
  ///     eps*^mu -> eps_T^mu = eps*^{mu nu} p_nu / M   (Eq.effectiveT)
  ///     m       -> m_T,
  /// so the whole spin-2 polarization tensor is consumed at the weak
  /// vertex through the single contracted object eps_T^mu, while the
  /// full tensor remains available to the subsequent strong decay.
  ///
  /// Note that eps_T is NOT normalised like a polarization vector; the
  /// extra 1/M in Eq.(effectiveT) is the note's convention and is part
  /// of the definition of the T_i form factors.
  ///
  /// EXPERIMENTAL STATUS: there is no standalone weak-current fit for
  /// D->K*_2(1430) analogous to {A_1(0),r_V,r_2}. The 2026 BESIII
  /// D^0 -> K pi e nu amplitude analysis FIXED the D-wave ratios
  ///     rVT = T_V(0)/T_1(0) = 1,   r2T = T_2(0)/T_1(0) = 1
  /// as model assumptions; they are NOT measurements. The defaults
  /// here reproduce that assumption and should be treated as such.
  //////////////////////////////////////////////////////////////////////////
  class VA_P_T : public VA_P_X_Base {
  private:
    FF_P_X * p_VT, * p_A0T, * p_A1T, * p_A2T;
    double   m_A1T_0, m_rVT, m_r2T;
    bool     m_useratios;
  public:
    VA_P_T(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name);
    ~VA_P_T() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
