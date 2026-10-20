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
  /// All THREE tags are registered against this class, and the topology
  /// is measured from how many of the current's legs are incoming rather
  /// than assumed:
  ///   VA_FF_0 : both legs incoming  - an annihilating pair, e.g. the
  ///             e+e- side of e+e- -> mu+mu-.
  ///   VA_F_F  : one in, one out     - a fermion line passing through,
  ///             e.g. tau -> nu_tau, or a t-channel lepton line.
  ///   VA_0_FF : both legs outgoing  - a produced pair, e.g. the lepton
  ///             pair of a semileptonic decay.
  /// A mismatched tag warns; it cannot change the answer, since one
  /// implementation serves all three. What DOES throw is a current whose
  /// two legs are not the same fermion line, which no later stage can
  /// detect and which silently computes a process that does not exist.
  ///
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
  public:
    /// Validate the tag against the topology actually supplied, and that
    /// the two legs really are one fermion line. Called via
    /// Current_Base::SetNIn, not from the constructor: which legs are
    /// incoming is context the constructor does not have.
    void CheckTopology() const override;
    /// Zero couplings (v = a = 0), e.g. a photon on a neutrino line. See
    /// Current_Base::Vanishes and the check in Exchange_ME.
    inline bool Vanishes() const override
    { return m_cR==Complex(0.,0.) && m_cL==Complex(0.,0.); }
    VA_0_FF(const ATOOLS::Flavour_Vector& flavs,
            const std::vector<int>& indices,
            const std::string& name);
    ~VA_0_FF() {}
    void SetModelParameters(struct GeneralModel model) override;
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti) override;
  };

  /// Aliases for the other two topologies. Same implementation; they
  /// exist only so that all three can carry honest, separately checkable
  /// tags. Trivial derived classes rather than further
  /// DEFINE_CURRENT_GETTERs on VA_0_FF, because that macro specialises
  /// Getter<...> on the class and so can only be used once per class.

  /// One leg in, one out: tau -> nu_tau, b -> c, a t-channel lepton line.
  class VA_F_F : public VA_0_FF {
  public:
    VA_F_F(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name) :
      VA_0_FF(flavs,indices,name) {}
  };

  /// Both legs in: an annihilating pair, e.g. the e+e- side of
  /// e+e- -> mu+mu-. Only reachable in a scattering process.
  class VA_FF_0 : public VA_0_FF {
  public:
    VA_FF_0(const ATOOLS::Flavour_Vector& flavs,
            const std::vector<int>& indices,
            const std::string& name) :
      VA_0_FF(flavs,indices,name) {}
  };
}

#endif
