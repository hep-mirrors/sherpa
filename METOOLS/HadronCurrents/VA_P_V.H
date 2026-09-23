#ifndef METOOLS_Current_Library_VA_P_V_H
#define METOOLS_Current_Library_VA_P_V_H

#include "METOOLS/HadronCurrents/VA_P_X_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// D -> V l nu, X with J^P = 1^-  (K*(892), rho, omega, phi).
  /// charm_semileptonic.tex Sec."D to V", Eqs.(PVvector),(PVaxial).
  ///
  /// The axial current is regrouped before implementation so that the
  /// q^2=0 constraint Eq.(A0constraint) is EXACT by construction
  /// rather than numerically hoped for (note item 2). Collecting the
  /// q^mu pieces of Eq.(PVaxial),
  ///
  ///   A^mu = (eps*.q) G(q^2)/q^2 q^mu
  ///          + (M+m) A_1(q^2) eps*^mu
  ///          - A_2(q^2) (eps*.q)/(M+m) P^mu,
  ///
  ///   G(q^2) = 2m A_0(q^2) - (M+m)A_1(q^2) + (M-m)A_2(q^2)
  ///          = 2m [ A_0(q^2) - A_3(q^2) ],
  ///
  /// using Delta/(M+m) = M-m. Since A_3(0)=A_0(0) is precisely
  /// Eq.(A0constraint), G(0)=0 identically and the 1/q^2 is regular.
  /// A_0(0) is therefore never read from the YAML: it is computed
  /// event by event from A_1(0), A_2(0) and the GENERATED m=sqrt(s_X),
  /// which keeps G(0)=0 exact even when the resonance mass runs.
  //////////////////////////////////////////////////////////////////////////
  class VA_P_V : public VA_P_X_Base {
  private:
    FF_P_X * p_V, * p_A0, * p_A1, * p_A2;
    /// Experimental normalisation convention: A_1(0), r_V=V(0)/A_1(0),
    /// r_2=A_2(0)/A_1(0). Applied to the form factors at setup.
    double   m_A1_0, m_rV, m_r2;
    bool     m_useratios;
  public:
    VA_P_V(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name);
    ~VA_P_V() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
