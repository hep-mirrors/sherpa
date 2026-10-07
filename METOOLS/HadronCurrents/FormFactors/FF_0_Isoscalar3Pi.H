#ifndef METOOLS_HadronCurrents_FF_0_Isoscalar3Pi_H
#define METOOLS_HadronCurrents_FF_0_Isoscalar3Pi_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"
#include <map>

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
  // THERE ARE TWO PARAMETERISATIONS, picked by the ff_model number: 100, the
  // default, is Eq. (10) below; 300 is a Khuri-Treiman form that replaces the
  // whole rho sector with a solved dispersion relation and runs only over
  // 0.454-1.060 GeV.  See the KHURI-TREIMAN MODE block further down.
  //
  // THE DEFAULT MODEL is Eq. (10) of hep-ph/0512180 (EPJC 47 (2006) 617):
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

    // CONTACT TERM -- a direct gamma* -> 3 pi coupling, added to the rho(770)
    // sum.  Off by default (m_contact_a = 0), so the fitted Eq. (10) is
    // unchanged unless it is asked for.
    //
    // CMD-2 (Budker preprint 2006-028) Eq. (39) writes the phi -> 3 pi Dalitz
    // amplitude as  |A_n a e^{i phi} + A_rhopi|^2,  with A_rhopi their sum of
    // three rho propagators over the pion pairings, unit coefficients, no
    // isospin factors -- our h770 -- and A_n a constant fixed by Eq. (43) so
    // that a is the contact amplitude as a FRACTION of the r.m.s. rho-pi
    // amplitude over the Dalitz plot.
    // Measured: a = 0.101 +- 0.044 +- 0.017, phi = -2.91 +- 0.14 +- 0.07
    // (CMD-2, 3.3 sigma); a = 0.104 +- 0.010 +- 0.020, phi = -2.47 +- 0.08
    // +- 0.08 (KLOE 2003).  The two agree on the size and differ on the phase.
    // Both phases in their Table 8 are already in CMD-2's convention: a
    // footnote records that KLOE's published A_dir = a e^{i phi} must be read
    // as a e^{-i phi}, so CMD-2 conjugated it.  Take both from Table 8 and do
    // not re-sign either.
    //
    // SIGN TRAP.  Their Eq. (42) propagator is
    //     1/D(Q^2) = 1/(Q^2/M^2 - 1 + i sqrt(Q^2) Gamma/M^2) = -BreitWigner,
    // so A_rhopi = -h770 and the constant that reproduces their amplitude
    // against OUR sum carries an extra pi.  m_contact_phase is stored in CMD-2's
    // convention and the pi is applied in Construct(), so the published number
    // can be typed straight into the run card.
    //
    // CAVEATS, all confirmed against the paper:
    //  * A_n = 7.52 is INHERITED from CMD-2's own preliminary analysis, not
    //    derived here, and the paper never states the sqrt(s) at which the
    //    Eq.(43) integral is evaluated (it moves 0.3% across their 1017-1021
    //    MeV band, so it does not matter).
    //  * Gamma_rho(Q^2) is NEVER DEFINED in the paper.  Their rho parameters
    //    are M = 775.8, Gamma = 150.3 (PDG 2004), common to all three charge
    //    states.  Reproducing A_n numerically gives 7.50 for
    //    Gamma(M^2/s)(p/p_M)^3, 7.62 for the GS-like (M/sqrt(s)) form, 7.74
    //    for bare p^3 and 7.31 for a constant width -- so 7.52 mildly favours
    //    the first, but the discrimination is weak.
    //  * h770 uses GounarisSakuraiM, different again, so a = 0.101 here is
    //    close to but not identical with their 10%.  The spread over line
    //    shapes is a few percent.  Pinning it properly needs the Dalitz
    //    comparison against BESIII_2018_I1703033.
    //
    // Their N_0 is a FREE parameter of the fit -- the efficiency-corrected
    // PRODUCED yield, never published -- not the 79577 observed events, and
    // the comparison runs through a 198x198 detector response matrix.  That
    // makes the shape-only reading stronger, not weaker: nothing in a or phi
    // constrains the rate.
    double  m_contact_a, m_contact_phase, m_contact_norm;

    // SHAPE-PRESERVING normalisation.  CMD-2 Eq.(39) reads
    //     N_k = (N_0/Z) Int_k dXdY |P+ x P-|^2 |A_n a e^{i phi} + A_rhopi|^2
    // with N_0 the OBSERVED event count and Z the integral over the whole
    // Dalitz plot -- so the fit is normalised to the data and a, phi are a
    // SHAPE measurement.  They say nothing about the rate.
    //
    // Added at amplitude level without this, the contact term moves the total
    // cross section by +17% at the phi, using parameters that carry no
    // information about it, which would break the fitted sigma(q^2).  So the
    // amplitude is rescaled by sqrt(I_off/I_on), both Dalitz integrals at the
    // same q^2: the contact term then redistributes events and leaves the rate
    // alone, which is what was measured.  Set Contact_preserve_norm = 0 to get
    // the unrescaled amplitude instead.
    bool    m_contact_preserve;
    //! Dalitz integral of |eps|^2 |amplitude|^2 at fixed q^2; contact on/off.
    double  DalitzIntegral(const double & q2,const bool & with_contact);
    //! sqrt(I_off/I_on), cached -- one integral per distinct q^2, so a
    //! fixed-energy run pays for it once.
    double  ContactNorm(const double & q2);
    std::map<double,double> m_norm_cache;
    Complex m_contact_c;
    int    m_i_plus, m_i_minus, m_i_zero; // which index is which pion

    // KHURI-TREIMAN MODE (FormFactor 300).  The whole rho sector -- all three
    // towers, and E and F with them -- is replaced by the dispersive basis of
    // KT3Pi_Tables.H,
    //     F = [A BW_om + B BW_phi + C BW_om' + D BW_om''](q^2)
    //         x Sum_ij [ F_a(s_ij; q^2) + r F_b(s_ij; q^2) ],
    // so the isoscalar bracket stays and plays the part of the subtraction
    // constant a(q^2).  Three consequences, all deliberate:
    //
    //  * E and F are dropped.  Inside the tabulated window, 0.454-1.060 GeV,
    //    s_ij never exceeds 0.40 M^2 of the rho(1450) or 0.29 M^2 of the
    //    rho(1700) even at the top, so those terms are smooth off-shell tails
    //    that the A--D bracket absorbs.  Above the window that stops being
    //    true -- F alone carries 50-80% of the cross section above 1.6 GeV --
    //    so the mode throws rather than run there.  Below 0.454 the amplitude
    //    is zero, with a warning: the solver stops converging as the Dalitz
    //    region closes towards the three-pion threshold at 0.419 GeV, so
    //    there is nothing to tabulate in the last 35 MeV, where the phase
    //    space has barely opened.
    //  * No contact term.  A second subtraction supplies the same freedom a
    //    contact term does, and having both would fit one degree of freedom
    //    twice, so asking for both throws.
    //  * No resonance registry.  The rho is not a parameter at all: it sits
    //    where the GKPY pi pi phase passes 90 degrees, 773.6 MeV, not at the
    //    particle database's mass.
    bool    m_kt, m_kt_warned;
    double  m_kt_r;       // second subtraction constant over the first [m_pi^-2]

    // Eq. (10) assembled from Sum_Term/Product_Term rather than written out as
    // algebra below.  Each entry names the invariant it is evaluated at, which
    // is what lets one product mix a function of q^2 with a function of a
    // two-pion sub-mass, and one sum hold the same rho at all three pairings.
    // Owns the terms; the propagators above stay owned by this class, since
    // the charged rho appears at two pairings and must not be deleted twice.
    Sum_Term * p_amp;      //!< with the contact term, if it is on
    Sum_Term * p_amp_nc;   //!< always without it, for the normalisation
    Sum_Term * MakeAmplitude(const bool & with_contact);
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
