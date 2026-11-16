#ifndef YFS_Main_Beta0_2_VPInsertion_H
#define YFS_Main_Beta0_2_VPInsertion_H
/*!
  \file Beta0_2_VPInsertion.H

  YFS: BETA0_2_VP_INSERTIONS vertex - the lepton-loop (e, mu, tau) insertion in the
  internal photon of the one-loop vertex of a space-like line, by the dispersive
  method. A photon carrying Pi(k^2) is a superposition of massive photons of mass^2 z
  with weight w(z) = -Im Delta alpha(z)/(pi z), so

    F_i^VP(q2) = int dz w(z) F_i(q2; photon mass sqrt z),

  with F_1 (on-shell renormalised) and F_2 of the one-loop vertex with a massive
  photon. The line's O(alpha^2) term X of Beta0_2_Vertex.H gains F1VP + F2VP R1
  (R1 of PauliFactors). Tabulated in yfs-helpers/vp-insertions (vertex_massive_photon.wl,
  tabulate_vp.wl, mkheader.py), with checks there: Schwinger F2(0) = 1/2, the IR slope
  nu A - 1, and the dispersion integral reproducing Delta alpha to 1e-20.

  Built for the PDG lepton masses (m_e = 0.51099895069e-3, m_mu = 0.1056583755, m_tau =
  1.77693 GeV; set them in PARTICLE_DATA): another line or loop mass gives 0, reported once.
*/
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace YFS {

  /*!
    F1VP + F2VP R1 of the line p -> pp (mass M), summed over the e, mu and tau loops,
    in (alpha/pi)^2 per unit loop charge squared; R1 from the other line k -> kp (mass
    m). 0 for a line mass other than m_e or m_mu, or q2 outside the table.
  */
  double VPInsertionVertex(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &pp, double M,
                           const ATOOLS::Vec4D &k, const ATOOLS::Vec4D &kp, double m);

  //! whether the table's lepton masses are the run's (checked once, reported)
  bool VPInsertionMassesMatch();

  /*!
    The lepton-loop insertion in the one-loop boxes of mu e -> mu e (direct + crossed),
    net of Delta alpha_lep(t) times the photonic boxes (which the Dyson factor already
    carries): N(s, t) in (alpha/pi)^2, |M|^2 relative to the Born, times Q_a Q_b
    (Beta0_2_VPBox_Table.H, yfs-helpers/vp-insertions/boxes). Tabulated at the nominal s of
    150 and 160 GeV muons on electrons at rest; another s uses the nearer one at the same
    t (a reduced Born after ISR: O(alpha) of an O(alpha^2) term). 0 unless fl is a mu e
    pair with the table's masses, or for t outside the table.
  */
  double VPInsertionBox(const ATOOLS::Flavour_Vector &fl, double s, double t);

}

#endif
