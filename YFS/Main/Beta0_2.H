#ifndef YFS_Main_Beta0_2_H
#define YFS_Main_Beta0_2_H

/*!
  \file Beta0_2.H

  The photonic two-loop beta_0^(2) (YFS: BETA0_2 for YFS.NLO, ACRAIC: BETA0_2
  for ACRAIC). With the IR-finite one-loop factor split as delta = delta_f + w
  (delta_f photonic, w the rest),

    delta^(2) = (1/2) delta_f^2 + delta_f w + NLL  [+ C L + K per vertex pair]

    delta_f = sum_{i<j} f_ij ,   f_ij = K_ij - i Im K_ij (c - 1/2) ,
    K_ij    = (alpha/2pi) (-Z_i Z_j) [C^A_nu(s_ij) - 1] ,
    NLL     = (alpha/pi) c_NLL sum_{i<j} [(Q_i^2 + Q_j^2) K_ij
                + (alpha/2pi) (-Z_i Z_j) (Q_i^2 - Q_j^2) ln(m_j/m_i)] ,
    c_NLL   = 3/32 - (3/4) zeta_2 + (3/2) zeta_3 ,

  Z_i = theta_i Q_i (+1 incoming, -1 outgoing), s_ij = (p_i +- p_j)^2 and
  c = ACRAIC: VIRT_HEL_IM_C; the c dependence cancels in |delta|^2 +
  2 Re delta^(2). The bracket is BETA0_2 two_loop, the rest of the massive
  form factor at leading power (hep-ph/0508127 eq. 21), L = ln(Q^2/m^2),
  C = c_NLL, K = 11/8 + 17 pi^2/32 - 9/4 zeta_3 - 2 pi^4/45 - pi^2 ln2/2.
  A vertex pair is a same-side f fbar pair, or one flavour in and out with
  Q^2 = -t > 100 m^2 when there is no same-side pair. The two-loop boxes are
  in Beta0_2_Boxes.H. Photonic only; alpha is YFS: BETA0_2_ALPHA.
*/

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "YFS/Main/YFS_Options.H"

#include <string>
#include <vector>

namespace YFS {

  //! c_NLL = 3/32 - (3/4) zeta_2 + (3/2) zeta_3
  double Beta02NLLCoefficient();
  //! C + K, in (alpha/pi)^2 Q^4
  double Beta02TwoLoopConstant();

  //! KKMC's KKbvir::CnuA, exact in the masses, s and t channel
  Complex Beta02CnuA(double svar, double m1, double m2);

  //! The settings of one column's beta_0^(2).
  struct Beta02_Config {
    beta02::code      mode  = beta02::none;
    beta02pairs::code pairs = beta02pairs::all;
    //! the two-loop vertex constant of a line through a space-like 2 -> 2
    beta02vertex::code crossedvertex = beta02vertex::exact;
    //! (1/2)(sum_d f_d)^2 (true) or sum_d (1/2) f_d^2 (false, as KKMC's CEEX2)
    bool   cross = true;
    //! the delta_f w product
    bool   weak  = true;
    //! (1/2) w^2: w exponentiated together with delta_f, so that
    //! delta^(2) = delta^2/2 + the logarithmic two-loop terms (BETA0_2_W_SQUARE)
    bool   wsquare = true;
    //! two_loop beyond O(alpha^2): the higher orders of exp(delta_f)(1 + NLL)
    //! (BETA0_2_HIGHER)
    beta02higher::code higher = beta02higher::none;
    //! the coupling of f_ij and the NLL
    double alpha = 0.;
    //! ACRAIC: VIRT_HEL_IM_C
    double imc   = 0.5;
    inline bool On() const { return mode != beta02::none; }
  };

  //! The analytic photonic pieces at one Born point.
  struct Photonic_Factors {
    bool ok = false;
    Complex deltaf = Complex(0., 0.);    //!< sum_pairs f_ij
    Complex nll    = Complex(0., 0.);    //!< sum_pairs NLL_ij
    Complex cst    = Complex(0., 0.);    //!< sum over vertex pairs of (alpha/pi)^2 Q^4 (C + K)
    std::vector<Complex> f, n, c;        //!< f_ij, NLL_ij and the constant per selected pair
    //! two-loop boxes per physical helicity and Born-weighted (Beta0_2_Boxes.H)
    std::vector<std::vector<int> > boxhel;
    std::vector<Complex> boxd2;
    double boxavg = 0.;
  };

  //! f_ij and NLL_ij of every selected charged pair of p (legs 0 .. nin-1
  //! incoming); ok is false if a charged leg is massless or nothing is charged
  Photonic_Factors PhotonicFactors(const ATOOLS::Flavour_Vector &fl,
                                   const ATOOLS::Vec4D_Vector &p,
                                   const Beta02_Config &cfg, size_t nin = 2);

  //! delta^(2) for the complex one-loop factor delta of one helicity
  Complex Beta02Delta(const Photonic_Factors &f, const Complex &delta,
                      const Beta02_Config &cfg);

  //! delta^(2) for a real delta, with every imaginary part dropped
  double Beta02DeltaReal(const Photonic_Factors &f, double delta,
                         const Beta02_Config &cfg);

  /*!
    KKMC's two-loop vertex factor from F1 = (alpha/pi) Q^2 cL/2:
      ll:       F1 + F1^2/2
      llnll:    ll + 2 c_NLL (alpha/pi) Q^2 F1
      exp:      exp(F1) (1 + 2 c_NLL (alpha/pi) Q^2 F1) - 1
      two_loop: llnll + (alpha/pi)^2 Q^4 (C + K)
  */
  Complex Beta02VertexFactor(const Complex &F1, double alpha, double Q2,
                             beta02::code mode);

  //! one line describing the configuration, for the run log
  std::string Beta02Describe(const Beta02_Config &cfg);

}

#endif
