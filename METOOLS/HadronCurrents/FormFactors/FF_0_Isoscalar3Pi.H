#ifndef METOOLS_HadronCurrents_FF_0_Isoscalar3Pi_H
#define METOOLS_HadronCurrents_FF_0_Isoscalar3Pi_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"

namespace METOOLS {

  ///////////////////////////////////////////////////////////////////////////
  //
  // The isoscalar electromagnetic form factor for  gamma* -> pi+ pi- pi0.
  //
  // NOT a tau current: a photon has no axial coupling here, so the whole
  // amplitude sits in the anomalous piece,
  //
  //     A^mu = eps^{mu al be ga} p_+al p_-be p_0ga  F(s; s_+-, s_+0, s_-0)
  //
  // and this class supplies F together with that epsilon structure, through
  // Current().
  //
  // THE MODEL is Eq. (10) of hep-ph/0512180 (EPJC 47 (2006) 617):
  //
  //   F = [A BW_om + B BW_phi + C BW_om' + D BW_om''] H_rho770
  //     +  E BW_phi                                   H_rho1450
  //     +  F BW_om''                                  H_rho1700
  //
  // each H being an all-plus sum over the three two-pion pairings.  That
  // all-plus sum is forced, not chosen: an I = 0 three-pion state is totally
  // antisymmetric in isospin and the momentum epsilon is totally
  // antisymmetric, so |A|^2 must come out symmetric under permuting the pions.
  //
  // The structure that matters is that this is a SUM OF PRODUCTS: E gives
  // rho(1450) to the phi alone and F gives rho(1700) to the omega'' alone.  A
  // factorised form -- one isoscalar function times one rho tower -- pairs
  // every isoscalar with every rho and cannot express that at any choice of
  // couplings.  Above ~1.9 GeV the two-pion sub-masses reach the
  // rho(1450)/rho(1700) poles while the A--D isoscalars are far off theirs, so
  // those two terms carry the tail, and a factorised model loses it.
  //
  // LINE SHAPES come from the shared registry, so rho(770) is whatever is
  // registered for it and the masses track the particle database -- with one
  // deliberate exception, below.
  //
  // THE ISOSCALAR POLES ARE MODEL PARAMETERS, NOT PARTICLE PROPERTIES.
  // Two separate reasons, both measured:
  //
  //  - omega(782) and phi(1020).  Their registered widths RUN, being built
  //    from a decay phase-space integral.  That is right near the pole and
  //    wrong far above it, where Gamma(s) keeps growing: by 1.6 GeV the
  //    running-width omega(782) propagator is a factor 3 smaller than a
  //    constant-width one and rotated in phase by 72 degrees.  Since A = 18.20
  //    dominates every other coupling by ~20x, that wrecks the isoscalar
  //    bracket across 1.3-1.6 GeV and no coupling can compensate -- a fit
  //    tried, and drove D to the wrong sign.
  //
  //  - omega(1420) and omega(1650).  hep-ph/0512180 fitted omega(1375) and
  //    omega(1631) with Gamma = 0.250/0.245; the particle database carries
  //    1.419/0.17 and 1.670/0.31.  Those are different objects: the paper's
  //    are parameters of its own amplitude.
  //
  // So all four are evaluated as FixedBreitWigner at model-level poles, which
  // default to the values fitted against Belle 2024.  The three rho stay on
  // the registry, where they belong.
  //
  // TUNING.  A--F are the paper's couplings refitted against Belle 2024 with
  // these line shapes; they cannot be the published numbers, which were fitted
  // against the paper's own inline propagators.  chi2/ndf = 1.80 over
  // 1.05-3.5 GeV, against 2.56 for a faithful transcription of the paper.
  //
  ///////////////////////////////////////////////////////////////////////////

  class FF_0_Isoscalar3Pi : public FormFactor_Base {
  private:
    // omega(782), phi(1020), omega(1420), omega(1650) -- all FixedBreitWigner
    // at model-level poles; see the note above on why they are not registry
    // line shapes.
    std::vector<Propagator_Base *> m_iso;
    double m_M[4], m_G[4];

    // rho(770), rho(1450), rho(1700), neutral and charged, from the registry.
    // The neutral shape takes the pi+ pi- pair and the charged one the two
    // mixed pairs.
    Propagator_Base * p_r770_0,  * p_r770_c;
    Propagator_Base * p_r1450_0, * p_r1450_c;
    Propagator_Base * p_r1700_0, * p_r1700_c;

    double m_c[6];                        // A .. F
    int    m_i_plus, m_i_minus, m_i_zero; // which index is which pion

    // Eq. (10) assembled from Sum_Term/Product_Term rather than written out as
    // algebra below.  Each entry names the invariant it is evaluated at, which
    // is what lets one product mix a function of q^2 with a function of a
    // two-pion sub-mass, and one sum hold the same rho at all three pairings.
    // Owns the terms; the propagators above stay owned by this class, since
    // the charged rho appears at two pairings and must not be deleted twice.
    Sum_Term * p_amp;
    void BuildAmplitude();

    void FixParameters(const FF_Parameters & params);
    void Construct();
  public:
    FF_0_Isoscalar3Pi(const FF_Parameters & params);
    ~FF_0_Isoscalar3Pi();

    //! The full hadronic current: F times the epsilon structure.
    ATOOLS::Vec4C Current(const ATOOLS::Vec4D_Vector & moms);
    //! The scalar F alone, for a caller that supplies its own epsilon.
    Complex operator()(const ATOOLS::Vec4D_Vector & moms);
  };
}

#endif
