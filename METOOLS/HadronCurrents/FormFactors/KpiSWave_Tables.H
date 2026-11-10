#ifndef METOOLS_HadronCurrents_FormFactors_KpiSWave_Tables_H
#define METOOLS_HadronCurrents_FormFactors_KpiSWave_Tables_H

#include <complex>

namespace METOOLS {
  typedef std::complex<double> Complex;
  // ------------------------------------------------------------------
  // The scalar K pi form factor fbar_0(s) = f_0(s)/f_+(0), I = 1/2, of
  //
  //   L. von Detten, F. Noel, C. Hanhart, M. Hoferichter, B. Kubis,
  //   Eur. Phys. J. C 81 (2021) 420, arXiv:2103.01966,
  //
  // for the four fits of their Table 3 to Belle tau -> K_S pi nu data.
  // Fit 3 is the reference: its source term is compatible with zero, so it
  // is the closest of the four to pure K pi final-state interaction, which
  // is what a K pi S-wave line shape inside a larger amplitude needs.  The
  // other three are its systematic spread.
  //
  // WHAT IT CONTAINS.  A two-channel (K pi, K eta') dispersive form factor:
  // an Omnes function built on the dispersive I = 1/2 S-wave phase of
  // Pelaez and Rodas (arXiv:1602.08404, CFD solution), dressed by a
  // two-resonance potential that brings in the K0*(1430) and K0*(1950).
  // So it carries the kappa = K0*(700) the way pi K scattering demands --
  // a pole deep in the complex plane, an Adler zero below threshold and
  // elastic unitarity -- which no Breit-Wigner can, plus the K0*(1430)
  // with K eta' inelasticity.  fbar_0(0) = 1 exactly; it is real and
  // analytic below threshold.
  //
  // WHY TABLES.  Each value needs an Omnes integral, two self-energy
  // integrals and a 2x2 solve.  No numerical table exists in the
  // literature, so these are generated from the paper's printed formulae
  // and parameters by a standalone script, make_KpiSWave_Tables.py, kept
  // with the fitting code rather than in Sherpa.  What it needs that the
  // papers do not state is recorded below.  The tables reproduce the
  // paper's own figure of fbar_0 to 1e-5 in |fbar_0| for sqrt(s) =
  // 0.64-1.0 GeV, 1e-4 to 1.3 GeV and 3e-3 up to m_tau, the residuals being
  // the rounding of the printed parameters.
  //
  // TRAPS for anyone rebuilding it from the papers:
  //  * von Detten et al. do not say which Pelaez-Rodas phase they used.
  //    The CFD one reproduces their figure; the UFD one misses by 3%.
  //  * The Adler zero of the input phase is misprinted in arXiv:1602.08404
  //    (Delta^2 for 4 Delta^2 under the root): as printed, fbar_0 comes out
  //    about 9% wrong.  The generator uses the corrected form.
  //  * The resonance masses in the potential (1.5708 and 2.133 GeV) are
  //    BARE parameters.  The physical K0*(1430) pole is 1408 - i180 MeV.
  //    Never replace them by PDG values.
  //
  // CONVENTIONS.  Isospin as in von Detten et al., K0bar pi-.  Masses those
  // of the input phase, m_pi = 0.13957, m_K = 0.496, m_eta = 0.547,
  // m_eta' = 0.95778 GeV, so the tabulated K pi threshold (0.63557 GeV) sits
  // 2.3 MeV above the physical K- pi+ one.  Between the two the table gives
  // its below-threshold value, which is real and finite.
  //
  // VALIDITY.  The form factor is constrained by tau data only below m_tau,
  // and above the K eta' threshold (sqrt(s) > 1.454 GeV) the four fits
  // disagree badly -- Fit 3 has a zero at 1.536 GeV, Fits 1 and 2 grow to
  // |fbar_0| ~ 8-20 above m_tau.  The paper calls that region "not well
  // constrained at all".  Treat the S-wave there as known only within the
  // spread of the four fits.
  //
  // THE GRID.  s from 0 to 3.3 GeV^2, in four segments split at the K pi,
  // K eta and K eta' thresholds.  Inside a segment s = lo + (hi - lo)
  // sin^2(theta) on a uniform theta grid, which turns the square-root cusp
  // at every threshold into a smooth function of theta; the lookup is a
  // four-point cubic in theta, separately in Re and Im.  81 + 601 + 601 +
  // 401 nodes per fit.
  //
  // ACCURACY of the lookup against the model, between nodes, at random s
  // and just above each threshold: 2e-5 of max |fbar_0| or better, for all
  // four fits.
  // ------------------------------------------------------------------

  extern const int    kKpiSWave_nfit;
  extern const int    kKpiSWave_nseg;
  extern const int    kKpiSWave_ntot;
  extern const double kKpiSWave_seg_lo[];
  extern const double kKpiSWave_seg_hi[];
  extern const int    kKpiSWave_seg_n[];
  extern const int    kKpiSWave_seg_offset[];
  extern const double kKpiSWave_re[];
  extern const double kKpiSWave_im[];

  //! fbar_0(s) of von Detten et al. Fit `fit` (1-4), s in GeV^2, 0 <= s <= 3.3.
  //  Outside that range, or for a fit that does not exist, it throws.
  Complex KpiSWave_FBar0(const int fit, const double s);
}

#endif
