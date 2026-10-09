#ifndef YFS_Main_Beta0_2_H
#define YFS_Main_Beta0_2_H

/*!
  \file Beta0_2.H

  The universal photonic two-loop beta_0^(2) at LL + NLL
  (YFS: BETA0_2 for YFS.NLO, ACRAIC: BETA0_2 for ACRAIC). Design and
  derivation: NOTES-two-loop-review-2026-10-01.md sec. 4 and 7 and
  NOTES-beta0-2-2026-10-01.md.

  With the IR-finite one-loop factor of the paper (eq. delta-hel) split as
    delta = delta_f + w ,
  delta_f the analytic photonic part and w the rest (weak loops, IR-finite
  boxes, VP), the two-loop IR-finite remainder on the Born-level amplitude is
    delta^(2) = (1/2) delta_f^2 + delta_f w + NLL ,
  with
    delta_f = sum_{i<j} f_ij ,
    f_ij    = K_ij - i Im K_ij (c - 1/2)
    K_ij    = (alpha/2pi) (-Z_i Z_j) [C^A_nu(s_ij) - 1]   (app-conventions eq. Kij)
    NLL     = (alpha/pi)^2 (c_NLL/2) sum_i Q_i^2 sum_{j!=i} (-Z_i Z_j)
                 [C^A_nu(s_ij) - 1 + ln(m_j/m_i)]
            = (alpha/pi) c_NLL sum_{i<j} [(Q_i^2+Q_j^2) K_ij
                 + (alpha/2pi)(-Z_i Z_j)(Q_i^2-Q_j^2) ln(m_j/m_i)] ,
    c_NLL   = 3/32 - (3/4) zeta_2 + (3/2) zeta_3 = 0.663...
  Z_i = theta_i Q_i (+1 incoming, -1 outgoing, as NLO_Base::LoopHelicityIR),
  s_ij = (p_i + p_j)^2 for two legs on the same side, (p_i - p_j)^2 otherwise,
  c = ACRAIC: VIRT_HEL_IM_C.

  - Re f_ij is the IR-finite photonic one-loop remainder in OUR B convention:
    the default subtraction (YFS: Dim_Reg true, BVV_full_eps) equals KKMC's
    SBvirt to 1e-11 (BVV_full_impl vs KKbvir::SBvirt), and vertex - SBvirt
    = (alpha/pi)(L - 1) = 2 Re K_ij for an s-channel pair (KKMC F1_1).
  - Im f_ij: our B_V subtracts i Im K_ij ln(mu^2 e^c/|s_ij|) (app-conventions
    eq. BV). The one-loop massive vertex has Im F1 = -(alpha/2)[ln(lambda^2/s)
    + 3/2] (analytic continuation of the spacelike vertex), so the photonic
    Im delta of a vertex pair is (alpha/2)(c - 3/2) = Im K_ij (1 - (c - 1/2)):
    KKMC's Im F1_1 at c = 1/2. f_ij shifts with c exactly as delta does,
    which is what makes the O(alpha^2) sum c-independent (paper sec. virtual;
    review sec. 1): |delta|^2 + 2 Re delta^(2) is invariant under a common
    shift delta, delta_f -> delta - i phi/2, delta_f - i phi/2.
  - LL: abelian exponentiation, (1/2) delta_f^2 = (1/2)(F_ini + F_fin)^2 for
    an s-channel 2 -> 2 process; KKMC's F1_2 - F1_1 carries
    (alpha/pi Q^2)^2 cL^2/8 = F1_1^2/2 per vertex (KKceex.cxx MakeF1ini).
  - NLL: KKMC's (alpha/pi Q^2)^2 cL (3/32 - 3/4 zeta_2 + 3/2 zeta_3) per
    vertex (KKceex.cxx MakeF1ini/MakeF1fin; hep-ph/0006359 eqs. 97-106 from
    Berends, van Neerven, Burgers NPB 297 (1988) 429). 16 c_NLL = 3/2 - 2 pi^2
    + 24 zeta_3, i.e. 4 c_NLL = 3/8 - pi^2/2 + 6 zeta_3, the C_F^2 coefficient
    of delta(1-x) in the two-loop non-singlet splitting function P_qq^(1)
    (Curci, Furmanski, Petronzio 1980) in (alpha/2pi)^2: the O(alpha^2 L)
    term of the massive form factor is the two-loop collinear anomalous
    dimension of each leg (massification, A_m = prod_i Z_i^(1/2) A_0;
    Penin hep-ph/0501120, Becher-Melnikov 0704.3582, Mitov-Moch
    hep-ph/0612149), (alpha/pi)^2 (c_NLL/2) per leg and per log. The per-leg
    sum above reduces to KKMC's c_NLL cL for one s-channel vertex; the
    distribution of the logs over the pairs is an O(alpha^2 L^0) choice.
  - delta_f w: the factorised product of the photonic and the non-photonic
    one-loop, ISR x weak, FSR x weak, IFI x photonic (YFS: BETA0_2_WEAK_PRODUCT).
    KKMC carries the vertex x EW-dressing part (EW-dressed Born times
    (1+F1ini)(1+F1fin)) but no box x F1; off reproduces KKMC.
  - O(alpha^2 L^0): BETA0_2 two_loop adds, per vertex pair, the rest of the
    exact photonic two-loop massive form factor at leading power (Penin
    hep-ph/0508127 eq. 21, Euclidean Q^2 = -q^2):
      f^(2) = (1/2) f^(1)^2 + C L + K ,   L = ln(Q^2/m^2) [= ln(s/m^2) - i pi],
    A vertex pair is a same-side f fbar pair (s channel), or one flavour in
    and out (Q^2 = -t, e.g. the electron line of mu e -> mu e) when no
    same-side vertex exists, no flavour is twice on one side and Q^2 > 100 m^2.
      C = c_NLL,  K = 11/8 + 17 pi^2/32 - 9/4 zeta_3 - 2 pi^4/45 - pi^2 ln2/2.
    The soft divergences exponentiate at one loop, so the YFS remainder is
    (1/2) r_1^2 + C L + K whatever the regulator and the split of B. KKMC's
    NLL is C (L - 1): two_loop = llnll + (alpha/pi)^2 Q^4 (C + K) per vertex.
    Other pairs (initial-final, final-state pairs of different currents)
    stay at LL + NLL: their constants are the two-loop boxes. No fermion
    loops (photonic only).

  The coupling: alpha of f_ij and of the NLL is YFS: BETA0_2_ALPHA (alpha0
  default, as KKMC). The one-loop delta of the provider is at the loop's
  coupling (G_mu in the paper cards): the YFS subtraction is divided by
  m_rescale_alpha = alpha(0)/alpha_model (NLO_Virtual.C), so its photonic part
  is (alpha_Gmu/alpha0) Re K_ij. With BETA0_2_WEAK_PRODUCT the difference
  ends up in delta_f w and delta^(2) = delta_f (delta - delta_f/2) depends on
  BETA0_2_ALPHA only at O((Delta r delta_f)^2) ~ 1e-6.
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
  //! C + K: the per-vertex O(alpha^2 L^0) of the exact photonic two-loop
  //! form factor beyond KKMC's C (L - 1), in (alpha/pi)^2 Q^4
  double Beta02TwoLoopConstant();

  //! KKMC's KKbvir::CnuA, exact in the masses, s and t channel (the same
  //! function as Acraic_Base::CnuA and NLO_LoopHelicity.C's LoopCnuA).
  Complex Beta02CnuA(double svar, double m1, double m2);

  //! The settings of one column's beta_0^(2).
  struct Beta02_Config {
    beta02::code      mode  = beta02::none;
    beta02pairs::code pairs = beta02pairs::all;
    //! (1/2)(sum_d f_d)^2 (true) or sum_d (1/2) f_d^2 (false, no cross
    //! products between pairs: KKMC's CEEX2 - CEEX1 step, which already has
    //! the ISR x FSR product at O(alpha))
    bool   cross = true;
    //! the delta_f w product
    bool   weak  = true;
    //! the coupling of f_ij and the NLL
    double alpha = 0.;
    //! ACRAIC: VIRT_HEL_IM_C, the phase convention of B_V
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
    //! BETA0_2 two_loop, e+e- -> l+l- only: the photon-exchange two-loop
    //! boxes per physical helicity (labels e-, e+, l-, l+) and Born-weighted
    //! (Beta0_2_Boxes.H); the per-helicity terms are added by the callers
    std::vector<std::vector<int> > boxhel;
    std::vector<Complex> boxd2;
    double boxavg = 0.;
  };

  /*!
    f_ij and NLL_ij of every selected charged pair of the point p (legs
    0 .. nin-1 incoming). ok is false if a charged leg is massless (the
    collinear logs need the mass) or nothing is charged.
  */
  Photonic_Factors PhotonicFactors(const ATOOLS::Flavour_Vector &fl,
                                   const ATOOLS::Vec4D_Vector &p,
                                   const Beta02_Config &cfg, size_t nin = 2);

  /*!
    delta^(2) for the complex one-loop factor delta = delta_f + w of one
    helicity (ACRAIC_Virtual: helicity).
  */
  Complex Beta02Delta(const Photonic_Factors &f, const Complex &delta,
                      const Beta02_Config &cfg);

  /*!
    delta^(2) with every imaginary part dropped: delta_f -> Re delta_f and a
    real delta (ACRAIC_Virtual: external, delta = v/2; YFS.NLO). The external
    delta has no Im part, so the m2 that goes with it must not have one
    either: with the complex delta_f it would keep (Im delta_f)^2, which
    depends on the phase convention c.
  */
  double Beta02DeltaReal(const Photonic_Factors &f, double delta,
                         const Beta02_Config &cfg);

  /*!
    KKMC's two-loop vertex factor from its one-loop F1 = (alpha/pi) Q^2 cL/2:
      ll:    F1 + F1^2/2
      llnll: F1 + F1^2/2 + 2 c_NLL (alpha/pi) Q^2 F1      (= KKMC F1_2)
      exp:   exp(F1) (1 + 2 c_NLL (alpha/pi) Q^2 F1) - 1
      two_loop: llnll + (alpha/pi)^2 Q^4 (C + K), the exact photonic vertex
    for the ACRAIC hand-coded virtual (ACRAIC_Virtual: acraic).
  */
  Complex Beta02VertexFactor(const Complex &F1, double alpha, double Q2,
                             beta02::code mode);

  //! one line describing the configuration, for the run log
  std::string Beta02Describe(const Beta02_Config &cfg);

}

#endif
