#ifndef METOOLS_Current_Library_VA_P_S_H
#define METOOLS_Current_Library_VA_P_S_H

#include "METOOLS/HadronCurrents/VA_P_X_Base.H"

namespace METOOLS {
  /// D -> S l nu, X with J^P = 0^+  (K*_0(1430), f0(980), f0(500),
  /// a0(980)). charm_semileptonic.tex Sec."D to S", Eq.(Scurrent):
  /// V^mu = 0 in this phase convention and
  ///   H^mu = +i { f+^S [P^mu - (Delta/q^2)q^mu]
  ///              + f0^S (Delta/q^2) q^mu }.
  /// The global factor i is a convention (it can be absorbed into a
  /// rephasing of the scalar state) but must then be used everywhere.
  class VA_P_S : public VA_P_X_Base {
  private:
    FF_P_X * p_fplus, * p_fzero;
    bool     m_tiefzero;
    /// Overall phase of the scalar current, exp(i*m_phase). Default
    /// pi/2, i.e. the "+i" of Eq.(Scurrent).
    double   m_phase;
  public:
    VA_P_S(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name);
    ~VA_P_S() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
