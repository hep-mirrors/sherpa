#ifndef EXTRAXS__NNLO__ee_aa_QED_TwoLoop_H
#define EXTRAXS__NNLO__ee_aa_QED_TwoLoop_H
/*!
  \file ee_aa_QED_TwoLoop.H

  Two-loop QED corrections to e+e- -> gamma gamma from Anastasiou, Glover,
  Tejeda-Yeomans, hep-ph/0201274 (AGT): spin-summed, massless electrons,
  conventional dimensional regularisation, finite remainders after Catani's
  I^(1), I^(2) subtraction, in the expansion
    A = 16 pi^2 alpha^2 [ T + (alpha/2pi) A^(1) + (alpha/2pi)^2 A^(2) ] ,
  T0 = 8 (t/u + u/t), x = -t/s, S = ln(s/mu^2), s-channel functions only.

  These are the scheme-independent CDR objects. The conversion to a massive,
  YFS-subtracted correction (massification, the YFS virtual form factor) is
  the user's: for ACRAIC/YFS BETA0_2 two_loop it is
  YFS::AddPhotonicTwoLoopDiphoton (YFS/Main/Beta0_2_Boxes.C).

  Polylogarithms of real argument are evaluated locally (ATOOLS has only
  DiLog). Checks and derivations: yfs-helpers/diphoton-two-loop (mpmath
  reference, Recola per-helicity, pole and scale checks).
*/
#include <complex>

namespace EXTRAXS {
  namespace ee_aa_QED_TwoLoop {

    //! real polylogarithms Li_n(w), n = 2, 3, 4, any real w <= 1
    double Li(int n, double w);

    // AGT App. A.1 / B.1, s-channel, raw (NOT divided by T0)
    double As (double x, double S = 0.);   //!< light-by-light, x sum_f Q_f^4 (eq. 4.7)
    double Bs (double x, double S = 0.);   //!< photonic ("abelian gluon")
    double E4s(double x, double S = 0.);   //!< x N_F'   (MSbar, massless loops)
    double F2s(double x, double S = 0.);   //!< x N_F'^2 (MSbar, massless loops)
    double G1s(double x, double S = 0.);   //!< one loop squared, photonic
    double X4s(double x, double S = 0.);   //!< x N_F'   (MSbar)
    double X5s(double x, double S = 0.);   //!< x N_F'^2 (MSbar)
    double T0 (double x);                  //!< 8 (t/u + u/t)

    //! Finite^{0x2}/T0, eq. (4.7)
    double Fin02(double x, double S = 0., double sumQ4 = 0., double NF = 0.);
    //! Finite^{1x1}/T0, eq. (5.4)
    double Fin11(double x, double S = 0., double NF = 0.);
    /*!
      The Catani-finite one-loop remainder <M0|M1>/T - I^(1), mu^2 = s, spin
      averaged, alpha/2pi (eqs. 3.6, 3.9-3.10, with the master integrals' overall
      sign fixed so that the poles cancel against I^(1)).
    */
    std::complex<double> R1(double x);

    //! the massless vertex form factor's one-loop Catani remainder (-4) and
    //! the eps^0 term f0 of ln F^(2) (C_F^2 part, alpha/4pi units, mu^2 = Q^2)
    double RF();
    double F0LogFinite();
  }
}

#endif
