#ifndef EXTRAXS__NNLO__ee_aa_QED_Poles_H
#define EXTRAXS__NNLO__ee_aa_QED_Poles_H
/*!
  \file ee_aa_QED_Poles.H

  Laurent series of the photonic one- and two-loop e+e- -> gamma gamma
  amplitudes (AGT, hep-ph/0201274), massless and massified, for the YFS
  double-virtual pole check (YFS NLO_Base::CheckVVPoles; derivation and
  references: yfs-helpers/diphoton-two-loop, massive_poles.py,
  yfs_pole_check.py).

  Conventions: D = 4 - 2 eps, mu^2 = s, coefficients in (alpha/2pi)^L, every
  quantity a ratio to the D-dimensional Born T(eps) (CDR; HV ratios to T0
  have the same poles). Normalisation:
    agt: MSbar, S_eps = (4pi)^eps e^{-eps gamma} per loop (AGT, MVV, Mitov-Moch)
    ol : (4pi)^eps / Gamma(1-eps) per loop (OpenLoops polenorm 0, Recola with
         delta_IR2 = pi^2/6, the finite part of Sherpa's BVV_full_eps);
         ol = agt x R^L, R = e^{-eps gamma} Gamma(1-eps) = 1 + zeta2/2 eps^2 + ...
  Massification: Mitov, Moch hep-ph/0612149 eqs. (36)-(39), C_F^2 part
  (photonic: C_A = n_f = 0), M_massive = Z M_massless with one Z for the two
  external electrons; Z^(1) through O(eps^2), Z^(2) through eps^0.
  Array index i holds the coefficient of eps^(i + lowest power), lowest power
  given in each comment.
*/
#include <complex>

namespace EXTRAXS {
  namespace ee_aa_QED_Poles {

    //! T(eps)/T0, eps^0 .. eps^2
    void TOverT0(double x, double t[3]);
    //! m(eps) = <M0|M1>/T(eps), complex, eps^-2 .. eps^2 (agt)
    void OneLoopSeries(double x, std::complex<double> m[5]);
    //! massless photonic two loop, per T(eps), eps^-4 .. eps^0 (agt):
    //! n02 = (P02 + Fin02)/T, n11 = (P11 + Fin11)/T; m1 = 2 Re m, eps^-2 .. eps^2
    void MasslessSeries(double x, double n02[5], double n11[5], double m1[5]);
    //! massified one loop m1m (eps^-1 .. eps^1) and two loop m2m (eps^-2 .. eps^0),
    //! L = ln(s/m^2) (mu^2 = s); ol = true: OpenLoops/Sherpa normalisation
    void MassiveSeries(double x, double L, double m1m[3], double m2m[3], bool ol = false);
  }
}

#endif
