#ifndef METOOLS_HadronCurrents_FormFactors_KT3Pi_Tables_H
#define METOOLS_HadronCurrents_FormFactors_KT3Pi_Tables_H

#include <complex>

namespace METOOLS {
  typedef std::complex<double> Complex;
  // ------------------------------------------------------------------
  // Basis functions for a dispersive (Khuri-Treiman) description of the
  // isoscalar gamma* -> pi+ pi- pi0 form factor: the alternative to the
  // Breit-Wigner sum of Eq. (10) of hep-ph/0512180 in FF_0_Isoscalar3Pi.
  //
  // WHAT IS TABULATED.  Khuri-Treiman writes the anomalous form factor as
  //     F(q^2; s,t,u) = F(s) + F(t) + F(u),
  // one function per two-pion channel, whose right-hand cut is fixed by
  // pi pi P-wave unitarity and whose left-hand cut is the feed-in from the
  // two crossed channels.  With two subtractions the solution is linear in
  // two constants,
  //     F(s) = a(q^2) [ F_a(s; q^2) + r F_b(s; q^2) ],
  // and these tables are F_a and F_b.  F_a is the solution seeded by the
  // Omnes function Omega(s), F_b the one seeded by s Omega(s), with
  // Omega(0) = 1.  a(q^2) and r are NOT in the tables; they belong to the
  // form factor.
  //
  // q^2 enters only through the crossed-channel feed-in, so the dependence
  // is weak: with the overall constant divided out, the shape in s moves by
  // about 11% from the bottom of the grid to the top.  The resonance
  // structure in q^2 -- omega, phi and above -- is carried by a(q^2), not
  // here.
  //
  // WHY TABLES.  Each q^2 slice is the converged solution of coupled
  // dispersion integrals, ten iterations and a few minutes per slice, none
  // of which can be done per event.  The lookup is also far cheaper than
  // the Breit-Wigner sum it stands in for -- a bilinear interpolation
  // against thirteen propagator evaluations -- and since q^2 is fixed per
  // event, an event only ever reads two adjacent slices.
  //
  // SOURCE.  The HISKP Khuri-Treiman solver
  // (github.com/HISKP-ph/khuri_treiman_solver), the code behind Stamen,
  // Isken, Kubis, Mikhasenko and Niehus, arXiv:2212.11767: its
  // V -> 3 pi (J^PC = 1--) solution with max_subs = 2, 10 iterations,
  // cutoff 1000 m_pi^2, and threshold/pseudo-threshold matching parameters
  // epsilon = 0.1, validity = 0.05, epsilon_tilde = 0.3,
  // validity_tilde = 0.2.
  //
  // PI PI PHASE.  GKPY, Garcia-Martin et al., arXiv:1102.2183, Eqs. (A7)
  // and (A8) with the constrained (CFD) parameters of their Table VII.
  // This is NOT the phase 2212.11767 used.  Theirs is the Bern phase, whose
  // chain of citations ends at a paper listed as in preparation and never
  // published.  GKPY is the published alternative that Niecknig, Kubis and
  // Schneider (arXiv:1203.2501) vary against it for exactly this purpose,
  // and they quote the difference as about 2% in Im F at the rho peak.
  //
  // VALIDATION.  At M = 7 m_pi the solver with this phase reproduces
  // |F|/|Omega| of Fig. 3 of 2212.11767 to 0.055 or better, against a
  // target digitised from that figure.  r = 0.03298 m_pi^-2 was fixed on
  // the BESIII omega -> 3 pi Dalitz parameter alpha.
  //
  // THE GRID.  78 solved values of M = sqrt(q^2)/m_pi from 3.2500 to
  // 7.5947, i.e. sqrt(q^2) = 0.454 to 1.060 GeV.  Spacing 0.0358
  // below M = 5.95 and 0.0717 above, plus two extra solves exactly at the
  // omega (5.6076) and the phi (7.3043), so at those two masses there is no
  // interpolation in q^2 at all.  In s, 200 uniform points from 4 m_pi^2
  // to (M_max - 1)^2 m_pi^2, the upper edge of the Dalitz plot at the top
  // mass.
  //
  // Every slice covers the WHOLE s range, including s above its own Dalitz
  // edge (M - 1)^2.  Those values are not extrapolated: the solver gives F
  // on the whole real axis, and interpolation in M near the edge of the
  // lower slice needs them.
  //
  // ACCURACY, measured against the solver output between nodes:
  //   in s, at a solved M:   5e-4 of |F|max pointwise, 1e-4 on the
  //                          Dalitz-integrated rate
  //   in M, at the spacing   2e-4 over most of the plot, 1.2e-3 at the
  //   above:                 rho, 6e-3 in the strip nearest the Dalitz
  //                          edge, 4e-4 on the rate
  // The edge is worst because each slice has a pseudo-threshold branch
  // point at s = (M - 1)^2, and interpolating in M at fixed s mixes one
  // slice's kink with its neighbour's.  The phase-space weight vanishes
  // there, which is why the rate barely notices.  All of it sits well
  // inside the ~2% phase-input systematic above.  The Dalitz parameters
  // alpha, beta, gamma computed from the tables agree with those from the
  // solver output to 1e-5 at the omega and 1e-3 at the phi.
  //
  // UNITS.  Everything is in the solver's units, m_pi = 0.13957 GeV: s is
  // in m_pi^2, F_b carries one power of s and so is in m_pi^2, and r is in
  // m_pi^-2.  KT3PiLookup takes q^2 and s in GeV^2 and converts, but what
  // it returns is in the solver's units.  Do not rescale F_b without
  // rescaling r.
  //
  // The numbers are the solver output resampled onto this grid.  Do not
  // edit them by hand; regenerate them from the solver.
  // ------------------------------------------------------------------

  const double kKT3Pi_mpi  = 0.13957;      // the solver's unit [GeV]
  const int    kKT3Pi_NM   = 78;         // q^2 slices
  const int    kKT3Pi_NS   = 200;        // s points per slice
  const double kKT3Pi_s_lo = 4.0;          // [m_pi^2]
  const double kKT3Pi_s_hi = 43.49006809;  // [m_pi^2], (M_max - 1)^2

  // M = sqrt(q^2)/m_pi of each slice, ascending and NOT uniformly spaced.
  extern const double kKT3Pi_M[78];

  // F_a and F_b, real and imaginary parts, slice-major: element
  // [j*kKT3Pi_NS + i] is slice j at s = s_lo + i (s_hi - s_lo)/(NS - 1).
  extern const double kKT3Pi_a_re[15600];
  extern const double kKT3Pi_a_im[15600];
  extern const double kKT3Pi_b_re[15600];
  extern const double kKT3Pi_b_im[15600];

  // Bilinear interpolation of one basis function: pass the _re and _im
  // arrays of F_a, or those of F_b.  q2 and s in GeV^2.
  //
  // Returns 0 when sqrt(q^2) lies outside the tabulated window rather than
  // clamping to the end slice.  Below 0.454 GeV that is a real gap, not a
  // kinematic impossibility -- three-pion production starts at 0.419 GeV
  // and an ISR run will ask for it -- so whether zero is acceptable there
  // is for the form factor to decide.
  //
  // s IS clamped onto the grid: every physical s lies inside it, so only
  // rounding at the Dalitz boundary can take s out.
  Complex KT3PiLookup(const double & q2,const double & s,
                      const double * re,const double * im);
}

#endif
