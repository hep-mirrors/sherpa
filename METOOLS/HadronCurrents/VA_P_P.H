#ifndef METOOLS_Current_Library_VA_P_P_H
#define METOOLS_Current_Library_VA_P_P_H

#include "METOOLS/HadronCurrents/VA_P_X_Base.H"

namespace METOOLS {
  /// D -> P l nu, X with J^P = 0^-  (K, pi, eta, eta').
  /// charm_semileptonic.tex Sec."D to P", Eq.(Pcurrent). A^mu = 0 by
  /// parity, so H^mu = V^mu = f+ [P^mu - (Delta/q^2) q^mu]
  ///                        + f0 (Delta/q^2) q^mu.
  class VA_P_P : public VA_P_X_Base {
  private:
    FF_P_X * p_fplus, * p_fzero;
    bool     m_tiefzero;   ///< impose f0(0)=f+(0) (on by default)
  public:
    VA_P_P(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name);
    ~VA_P_P() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
