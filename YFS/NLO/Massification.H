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
    YFS: VIRTUAL_MASSIFICATION (off, electrons, leptons, light_fermions,
    light_fermions_charm; default off). The one-loop provider is called with
    the selected external species (e; e, mu; e, mu, u, d, s; or those and c)
    massless, at the massless image of the event's point, and the small-mass
    amplitude is rebuilt with the per-leg massification factor (1811.06461)
    at the |M|^2 level

      delta_massive = delta_massless
                    + sum_legs Q^2 [L^2 + L + 4 + zeta2 + c_provider],

    delta = 2 Re(A0^* A1)/|A0|^2 in units of alpha/4pi, L = ln(mu_IR^2/m^2),
    c_provider the provider's 1/eps normalisation (Virtual_ME2_Base::
    MasslessLegConstant). Exact up to O(m^2/Q^2) for every invariant Q^2.

    Fails at a pair threshold (excluded at setup, 4 m^2/s), for |t| < ~100 m^2
    and for a photon inside the dead cone of a massless lepton in the
    real-virtual (both counted per run and reported at the end).
    The tau is never massified. Closed loops of the massless species are
    massless too, a power correction m^2/|q^2|.
  */
  class Massification {
  public:
    static bool Enabled();
    static virtmassif::code Mode();
    //! flavs: the loop process, incoming first. Registers its massified
    //! species as massless; call before the provider's GetME2.
    Massification(const ATOOLS::Flavour_Vector &flavs, const std::string &what);
    //! after GetME2: the provider's normalisation constant (throws if it keeps the masses)
    void Attach(const PHASIC::Virtual_ME2_Base *me);
    /*!
      The massless image of p: in the incoming pair's rest frame, directions
      are kept and the outgoing three-momenta rescaled by a common factor so
      the energies (provider masses) add up to sqrt(s); then boosted back.
    */
    ATOOLS::Vec4D_Vector MasslessImage(const ATOOLS::Vec4D_Vector &p) const;
    //! sum_legs Q^2 z_i, in units of alpha/2pi (the providers' ME_Finite)
    double Z(const double muir2) const;
    //! the 1/eps coefficient, sum_legs Q^2 (2L + 1)/2 in alpha/2pi
    double ZPole1(const double muir2) const;
    //! the 1/eps^2 coefficient, sum_legs Q^2 in alpha/2pi (cancels the provider's double pole)
    double ZPole2() const;
    //! checks that the provider's 1/eps^2 (ME_E2 units) + ZPole2 vanishes; counts failures
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
    //! the provider takes the small-mass limit itself (Recola): no Z factor or double-pole check
    bool m_native = false;
  };

}

#endif
