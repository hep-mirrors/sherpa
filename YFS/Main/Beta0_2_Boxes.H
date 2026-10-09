#ifndef YFS_Main_Beta0_2_Boxes_H
#define YFS_Main_Beta0_2_Boxes_H
/*!
  \file Beta0_2_Boxes.H

  The non-factorisable photonic two-loop beta_0^(2) of e+e- -> l+l-
  (l = mu, tau), part of BETA0_2 two_loop: the two-loop QED amplitude
  beyond the two vertex form factors, for photon exchange.

  Massification (A_m = prod_legs Z_i^(1/2) A_0, photonic: no soft function)
  applied to the amplitude and to each vertex form factor gives, in log form,

    ln A_m - B/2 = sum_vertex [ln F_m - B/2] + [ln A_0 - ln F_0 - ln F_0 - B_IF/2] .

  The first term is the exact massive vertex (Beta0_2.C, C L + K). The second
  has no two-loop poles (photonic, Nf = 0: no two-loop cusp, beta_0 = 0; the
  two-loop collinear poles are per leg and cancel), and B_IF is one-loop, so
  its two-loop part is the massless, IR-scheme independent

    N(h) = [R2/R0 - (R1/R0)^2/2] - 2 [ln R^F]^(2)          (alpha/2pi)^2 ,

  R the IR-subtracted helicity amplitudes of arXiv:2411.07063 (mu^2 = s) and
  [ln R^F]^(2) = 4 (c_a - i pi C) the massless form factor in the same scheme
  (Penin hep-ph/0508127 eq. 20). N is tabulated in x = -t/s
  (Beta0_2_Boxes_Table.H) and enters delta^(2)(h) as

    (alpha/2pi)^2 N(h, x) A_gamma(h)/A_0(h) ,

  the photon-exchange share of the helicity amplitude (massless chiral
  couplings, as Dyson_VP). Not included: the two-loop boxes with a Z
  exchanged, and anything for other final states (the source has equal
  charges, so quarks cannot be separated).
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
    The per-helicity (physical labels, legs e-, e+, l-, l+) and the
    Born-weighted two-loop box terms of the point p, added to f. No-op unless
    fl is e- e+ -> l- l+ with l = mu or tau.
  */
  void AddPhotonicTwoLoopBoxes(Photonic_Factors &f,
                               const ATOOLS::Flavour_Vector &fl,
                               const ATOOLS::Vec4D_Vector &p, double alpha);

  //! the box term of the physical helicity labels l (0 if none)
  Complex Beta02BoxFor(const Photonic_Factors &f, const std::vector<int> &l);

  /*!
    e+e- -> gamma gamma only (no-op otherwise): the exact photonic two-loop
    hard function of hep-ph/0201274, read from the double-virtual provider
    vv (VV_Generator: Internal, EXTRA_XS/NNLO/ee_aa_QED_VV.H; 0 for any
    other provider or none), evaluated at p, massified, as a
    helicity-independent real correction to the BETA0_2 two_loop delta^(2)
    of the e+e- pair:
      f.cst += (1/2) (alpha/2pi)^2 Delta(x),   x = -t/s.
    Assumes BETA0_2_WEAK_PRODUCT true. Returns the added term; provided:
    vv is the gamma gamma provider; skipped: m_e^2/min(|t|,|u|) > powercut
    (beyond leading power), nothing added.
  */
  double AddPhotonicTwoLoopDiphoton(Photonic_Factors &f,
                                    PHASIC::Virtual_ME2_Base *vv,
                                    const ATOOLS::Vec4D_Vector &p, double alpha,
                                    double me, double powercut,
                                    bool &provided, bool &skipped);

}

#endif
