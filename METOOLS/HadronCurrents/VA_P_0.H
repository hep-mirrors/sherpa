#ifndef METOOLS_Current_Library_VA_P_0_H
#define METOOLS_Current_Library_VA_P_0_H

#include "METOOLS/HadronCurrents/Current_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// Pseudoscalar decay-constant current for a PURELY LEPTONIC decay,
  ///
  ///     <0| qbar gamma^mu (1-gamma_5) c |P(p)> = i f_P V_qq' p^mu,
  ///
  /// i.e. D+ -> mu+ nu, D_s+ -> tau+ nu, pi+ -> mu+ nu.
  ///
  /// NAMING. The legacy tag was VA_0_P, but by this directory's
  /// <in>_<out> convention that reads as "nothing in, one pseudoscalar
  /// out", which is backwards: the pseudoscalar is the DECAYING
  /// particle and annihilates into the current. Every use is P -> 0,
  /// so the class is VA_P_0. VA_0_P stays registered as an alias, with
  /// the same treatment as VA_F_F / VA_0_FF: the topology is checked
  /// against the tag rather than assumed, so a mislabel is visible.
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
    void SetModelParameters(struct GeneralModel model);
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
  };

  /// Deprecated alias for the legacy tag. Same class.
  class VA_0_P : public VA_P_0 {
  public:
    VA_0_P(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name) :
      VA_P_0(flavs,indices,name) {}
  };
}

#endif
