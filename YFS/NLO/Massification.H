#ifndef YFS_NLO_Massification_H
#define YFS_NLO_Massification_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "YFS/Main/YFS_Options.H"
#include <string>
#include <vector>

namespace PHASIC { class Virtual_ME2_Base; }

namespace YFS {

  /*!
    YFS: VIRTUAL_MASSIFICATION (off, electrons, leptons, light_fermions;
    default off). The one-loop provider is called with the external electrons
    (electrons), electrons and muons (leptons), or electrons, muons and the
    u, d, s quarks (light_fermions) massless, at the massless
    image of the event's point, and the small-mass amplitude is rebuilt with
    the per-leg massification factor (Penin; Mitov, Moch; Becher, Melnikov;
    Engel, Gnendiger, Signer, Ulrich 1811.06461), at the |M|^2 level

      delta_massive = delta_massless
                    + sum_legs Q^2 [L^2 + L + 4 + zeta2 + c_provider],

    delta = 2 Re(A0^* A1)/|A0|^2 in units of alpha/4pi, L = ln(mu_IR^2/m^2),
    c_provider the provider's 1/eps normalisation (Virtual_ME2_Base::
    MasslessLegConstant). Exact up to O(m^2/Q^2) for every invariant Q^2.

    Validated point by point (2026-10-06, full EW, G_mu, OpenLoops and
    Recola 2): residual <= 6e-4 alpha/4pi at the Z pole and above for
    e+e- -> mu mu (every helicity); it fails at a pair threshold (Coulomb,
    2 pi^2/beta), for |t| < ~100 m^2 (forward Bhabha) and for a photon
    inside the dead cone of a massless lepton in the real-virtual (3 - 10
    alpha/4pi there). The first is excluded at setup (4 m^2/s), the other two
    are counted per run and reported at the end.

    tau is never massified: m_tau^2/s is 4e-4 at the Z pole, where the
    residual would be ~0.5 alpha/4pi. Closed loops of the massless species
    are massless too (the provider's masses are global): a power correction
    m^2/|q^2| of the loop's virtuality.
  */
  class Massification {
  public:
    static bool Enabled();
    static virtmassif::code Mode();
    //! flavs: the loop process, incoming first. Registers its external
    //! massified species (see Mode) as massless loop species; call BEFORE
    //! the provider's GetME2.
    Massification(const ATOOLS::Flavour_Vector &flavs, const std::string &what);
    //! after GetME2: the provider's normalisation constant; throws if the
    //! provider keeps the masses
    void Attach(const PHASIC::Virtual_ME2_Base *me);
    /*!
      The massless image of p: in the rest frame of the incoming pair, the
      incoming momenta keep their directions with the provider's masses, the
      outgoing three-momenta are rescaled by one common factor so that the
      energies (provider masses: 0 for the massified legs, unchanged
      otherwise) add up to sqrt(s); then boosted back. Momentum conservation
      is exact, every direction in that frame is kept (inverse RAMBO).
    */
    ATOOLS::Vec4D_Vector MasslessImage(const ATOOLS::Vec4D_Vector &p) const;
    //! sum_legs Q^2 z_i, in units of alpha/2pi (the providers' ME_Finite)
    double Z(const double muir2) const;
    //! the 1/eps coefficient of the same factor, sum_legs Q^2 (2L + 1)/2 in
    //! alpha/2pi: massive pole = the provider's (massless) pole + this
    double ZPole1(const double muir2) const;
    //! the 1/eps^2 coefficient, sum_legs Q^2 in alpha/2pi: it cancels the
    //! provider's massless double pole leg by leg (massive pole = 0)
    double ZPole2() const;
    //! every evaluation: the provider's 1/eps^2 (ME_E2 units) + ZPole2 must
    //! vanish; counts and reports the events where it does not
    void CheckDoublePole(const double providere2);
    //! per-evaluation diagnostics on the massive point p
    void Count(const ATOOLS::Vec4D_Vector &p);
    //! end of run: the counters summed over the MPI ranks (NLO_Base)
    void Report();
  private:
    std::string m_what;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<size_t> m_legs, m_photons;
    std::vector<double> m_provmass;
    double m_c;
    long m_n, m_nforward, m_ndeadcone;
    long m_ne2bad = 0;
    double m_e2worst = 0.;
  };

}

#endif
