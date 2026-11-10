#ifndef METOOLS_Current_Library_VA_P_0_H
#define METOOLS_Current_Library_VA_P_0_H

#include "METOOLS/HadronCurrents/Current_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// Pseudoscalar decay-constant current, in either direction:
  ///
  ///   VA_P_0, P annihilated:  <0| qbar gamma^mu (1-gamma_5) q' |P(p)>
  ///                             = i f_P V_qq' p^mu
  ///           D+ -> mu+ nu, D_s+ -> tau+ nu, pi+ -> mu+ nu.
  ///   VA_0_P, P created:      <P(p)| qbar gamma^mu (1-gamma_5) q' |0>
  ///                             = -i f_P V_qq' p^mu
  ///           tau- -> pi- nu, tau- -> K- nu.
  ///
  /// NAMING. By this directory's <in>_<out> convention VA_P_0 is "one
  /// pseudoscalar in, nothing out" and VA_0_P "nothing in, one
  /// pseudoscalar out". The two are the same matrix element up to complex
  /// conjugation, so they are one class. As with VA_F_F / VA_0_FF, the tag
  /// is checked against the topology the current is handed, not assumed, so
  /// a mislabel stops the run.
  ///
  /// The decay constants are the legacy values, kept so that existing
  /// tables reproduce their previous numbers exactly. They are all
  /// overridable with fP, and V_qq' with Vxx.
  //////////////////////////////////////////////////////////////////////////
  class VA_P_0 : public Current_Base {
  private:
    double m_fP, m_Vxx;
  public:
    VA_P_0(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name) :
      Current_Base(flavs,indices,name), m_fP(1.), m_Vxx(1.) {}
    ~VA_P_0() {}
    void CheckTopology() const override;
    void SetModelParameters(struct GeneralModel model);
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
  };

  /// The created pseudoscalar, tau- -> P- nu. Same class; the tag only
  /// changes which topology CheckTopology accepts.
  class VA_0_P : public VA_P_0 {
  public:
    VA_0_P(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name) :
      VA_P_0(flavs,indices,name) {}
  };
}

#endif
