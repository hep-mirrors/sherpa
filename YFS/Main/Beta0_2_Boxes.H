#ifndef YFS_Main_Beta0_2_Boxes_H
#define YFS_Main_Beta0_2_Boxes_H
/*!
  \file Beta0_2_Boxes.H

  The non-factorisable photonic two-loop beta_0^(2) of e+e- -> l+l-
  (l = mu, tau), part of BETA0_2 two_loop: the two-loop QED amplitude
  beyond the two vertex form factors, for photon exchange.

  Massification applied to the amplitude and to each vertex form factor leaves,
  beyond the exact massive vertices (Beta0_2.C), a massless, IR-scheme
  independent two-loop remainder

    N(h) = [R2/R0 - (R1/R0)^2/2] - 2 [ln R^F]^(2)          (alpha/2pi)^2 ,

  R the IR-subtracted helicity amplitudes of arXiv:2411.07063 (mu^2 = s) and
  [ln R^F]^(2) the massless form factor in the same scheme (hep-ph/0508127).
  N is tabulated in x = -t/s (Beta0_2_Boxes_Table.H) and enters delta^(2)(h) as

    (alpha/2pi)^2 N(h, x) A_gamma(h)/A_0(h) ,

  the photon-exchange share of the helicity amplitude (massless chiral
  couplings). Not included: boxes with a Z exchanged, other final states.

  Bhabha uses the same N per helicity class (Beta0_2_Boxes_Bhabha_Table.H):
  s-channel only, s + t, and t-channel only. A_gamma is the s + t photon
  amplitude of the class. Applied at leading power only, min(|t|, |u|) > 100 m_e^2.
*/
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include <vector>

namespace PHASIC { class Virtual_ME2_Base; }

namespace YFS {

  struct Photonic_Factors;

  //! N(h, x) in (alpha/2pi)^2: same = e- and l- of the same helicity
  Complex TwoLoopPhotonicN(bool same, double x);

  /*!
    The per-helicity (physical labels, legs e-, e+, l-, l+) and Born-weighted
    two-loop box terms at p, added to f. No-op unless fl is e- e+ -> l- l+
    with l = e, mu or tau.
  */
  void AddPhotonicTwoLoopBoxes(Photonic_Factors &f,
                               const ATOOLS::Flavour_Vector &fl,
                               const ATOOLS::Vec4D_Vector &p, double alpha);

  //! the box term of the physical helicity labels l (0 if none)
  Complex Beta02BoxFor(const Photonic_Factors &f, const std::vector<int> &l);

  /*!
    e+e- -> gamma gamma only (no-op otherwise): the photonic two-loop hard
    function of hep-ph/0201274 from the double-virtual provider vv
    (VV_Generator: Internal), massified, as a helicity-independent correction
    to the BETA0_2 two_loop delta^(2):
      f.cst += (1/2) (alpha/2pi)^2 Delta(x),   x = -t/s.
    Assumes BETA0_2_WEAK_PRODUCT true. Returns the added term; provided: vv is
    the gamma gamma provider; skipped: m_e^2/min(|t|,|u|) > powercut, nothing added.
  */
  double AddPhotonicTwoLoopDiphoton(Photonic_Factors &f,
                                    PHASIC::Virtual_ME2_Base *vv,
                                    const ATOOLS::Vec4D_Vector &p, double alpha,
                                    double me, double powercut,
                                    bool &provided, bool &skipped);

}

#endif
