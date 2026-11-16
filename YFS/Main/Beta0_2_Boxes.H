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
    l e -> l e (l = mu or tau, e massless at leading power; BETA0_2_MUE_BOXES
    photonic): the photonic two-loop remainder beyond the vertex form factors,

      D = (1/2) M2/B - (M1/(2B))^2     on the three mixed charge classes
                                        (Q_e^3 Q_l, Q_e^2 Q_l^2, Q_e Q_l^3),

    B, M1 = 2 Re<M0|M1>, M2 = 2 Re<M0|M2> + <M1|M1> summed over spins, in CDR.
    D is IR finite and unchanged by any helicity-independent shift of the one
    loop (IR scheme, electron mass), so it needs neither subtraction nor
    massification. Source: the two-loop matrix element of arXiv:2212.06481
    (ancillary matrix-element.m, parsed as data, GPLs with handyG), tabulated
    in Beta0_2_Boxes_MuE_Table.H; the vertex classes stay with Beta0_2_Vertex.
    With the Born-weighted delta^(2) the missing O(alpha^2) is 2 D, so
    f.muebox = (alpha/pi)^2 D is added to f.boxavg. This needs BETA0_2_W_SQUARE
    (the full delta^2/2), and two_loop then leaves out the NLL terms of the
    four mu-e pairs, which are mixed-class and inside D.

    Applied only where the table holds (otherwise f.muebox = 0 and the status
    says why):
    - the Born is l e -> l e, elastic, l = mu or tau, e = e+ or e-, nothing else;
    - s/M^2 within the table's nodes (3 .. 16.5; MUonE 150/160 GeV are nodes);
    - -t > 100 m_e^2 (leading power in m_e^2/|t|, as for the Bhabha boxes),
      which also keeps t/tmin above the table's lowest node.
    Above the highest node (t -> tmin, theta_e -> 0) D has saturated and the
    edge value is used.
  */
  enum class MuEBox { not_mue, applied, outside_s, below_power };
  MuEBox AddMuEPhotonicTwoLoop(Photonic_Factors &f, const ATOOLS::Flavour_Vector &fl,
                               const ATOOLS::Vec4D_Vector &p, double alpha);

  //! D in (alpha/pi)^2 at s = s/M^2 > 1, t = t/M^2 < 0 for the charge product qQ = Q_e Q_l
  //! (table interpolation, clamped to the grid)
  double MuEPhotonicTwoLoopD(double s, double t, double qQ);

  /*!
    The l e -> l e term for a per-helicity delta^(2): |1 + delta_h + d2_h|^2
    already carries the helicity spread of the one loop, which D contains, so
    every helicity gets f.muebox - Var_h(Re delta_h), Born-weighted with
    |a0_h|^2 (the muon-vertex part of the spread, ~5e-7 for MUonE, is left in).
    0 if f.muebox is 0 or the factors are unusable.
  */
  double MuEBoxHelicityShare(const Photonic_Factors &f, const std::vector<Complex> &a0,
                             const std::vector<Complex> &delta);

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
