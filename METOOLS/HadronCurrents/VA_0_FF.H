#ifndef METOOLS_Current_Library_VA_0_FF_H
#define METOOLS_Current_Library_VA_0_FF_H

#include "METOOLS/HadronCurrents/Current_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// Elementary fermion-antifermion V-A current
  ///
  ///     J^mu = ubar(p_0) gamma^mu [ v - a gamma_5 ] u(p_1).
  ///
  /// NAMING. This directory's convention is <in>_<out>: VA_0_PP is
  /// "nothing in, two pseudoscalars out", VA_P_V is "pseudoscalar in,
  /// vector out". The legacy tag VA_F_F therefore reads as "one
  /// fermion in, one fermion out" - which is right for a tau -> nu_tau
  /// transition current, and wrong for the lepton pair of a
  /// semileptonic decay, where BOTH legs are outgoing. The class is
  /// named VA_0_FF after the latter, which is the topology the charm
  /// table uses throughout.
  ///
  /// The distinction is not cosmetic. The two topologies need
  /// different spinor types (an incoming u versus an outgoing v), and
  /// the existing decay tables already exercise both through this one
  /// class: the tau table has 30 entries of the form Indices: [N, 0] -
  /// index 0 being the decaying particle, i.e. genuinely F -> F - and
  /// 2 entries Indices: [3, 2] with both legs outgoing. XYZFunc
  /// resolves the spinor types itself, so a single implementation is
  /// correct for both; what was missing was any check that the tag
  /// matches the topology actually supplied.
  ///
  /// Both tags are therefore registered against this class and the
  /// topology is now VALIDATED rather than assumed:
  ///   VA_0_FF : neither index may be the decaying particle. Throws.
  ///   VA_F_F  : exactly one index should be the decaying particle.
  ///             Warns if not, rather than throwing, so that existing
  ///             tables keep running - but those entries are
  ///             mislabelled and should be switched to VA_0_FF.
  /// Index order is unchanged: 0 = the BARRED spinor, 1 = the unbarred
  /// one. Parameters: v (default 1), a (default -1).
  ///
  /// Ported from HADRONS++/Current_Library/VA_F_F.{H,C}. Two latent
  /// bugs were fixed rather than transcribed: p_ff was never
  /// initialised although the destructor unconditionally deleted it
  /// (the dead FF_Base/Simple_FF scaffolding is gone, which removes it
  /// by construction), and an unrecognised V_A_FORM_FACTOR called
  /// Abort() instead of throwing.
  //////////////////////////////////////////////////////////////////////////
  class VA_0_FF : public Current_Base {
  private:
    Complex m_cR, m_cL;
    /// True if one of the two legs is the decaying particle
    /// (external index 0), i.e. the F -> F topology.
    bool    m_hasincoming;
    void CheckTopology() const;
  public:
    VA_0_FF(const ATOOLS::Flavour_Vector& flavs,
            const std::vector<int>& indices,
            const std::string& name);
    ~VA_0_FF() {}
    void SetModelParameters(struct GeneralModel model);
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
  };

  /// Alias for the F -> F topology (tau -> nu_tau, b -> c, ...). Same
  /// implementation; it exists only so that the two topologies can
  /// carry honest, separately checkable tags. A trivial derived class
  /// rather than a second DEFINE_CURRENT_GETTER on VA_0_FF, because
  /// that macro specialises Getter<...> on the class and so can only
  /// be used once per class.
  class VA_F_F : public VA_0_FF {
  public:
    VA_F_F(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name) :
      VA_0_FF(flavs,indices,name) {}
  };
}

#endif
