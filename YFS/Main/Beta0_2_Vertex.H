#ifndef YFS_Main_Beta0_2_Vertex_H
#define YFS_Main_Beta0_2_Vertex_H
/*!
  \file Beta0_2_Vertex.H
  The exact photonic two-loop vertex of a line that passes through a 2 -> 2
  space-like scattering unchanged (the muon and the electron line of
  mu e -> mu e), used by BETA0_2 two_loop with YFS: BETA0_2_CROSSED_VERTEX
  exact (Beta0_2.C).

  A line of mass m and momentum transfer q2 = -t/m^2 > 0 has the on-shell
  vertex F1 gamma^a + F2 i sigma^{ab} q_b/(2m). With rho = F2/F1, which is IR
  finite, the spin-summed |M|^2 is |F1|^2 (T0 + 2 rho T1 + rho^2 T2). After the
  YFS subtraction its O(alpha^2) part is
     2 delta^2 + 2 X,
     X = [ln F1]^(2) + rho^(2) R1 + rho^(1)^2 (R2 - 2 R1^2)/2 ,
  where delta is the one-loop coefficient (it already contains rho^(1) R1) and
  R_i = T_i/T0 (PauliFactors). X, in units of (alpha/pi)^2 Q^4, replaces the
  constant C L + K of hep-ph/0508127, which is its q2 -> infinity limit.

  The remainders are tabulated in ln q2 (Beta0_2_Vertex_Table.H, from the
  photonic part of arXiv:1712.09889) for 1e-6 < q2 < 1e8. Below, [ln F1]^(2) is
  continued linearly in q2 and rho at its q2 -> 0 value; above, [ln F1]^(2)
  follows C L + K and rho is set to its last value.
*/
#include "ATOOLS/Math/Vector.H"

namespace YFS {

  //! [ln F1]^(2) = F1^(2) - (F1^(1))^2/2, in (alpha/pi)^2, q2 = -t/m^2 > 0
  double TwoLoopVertexLnF1(double q2);
  //! rho^(1) = F2^(1), in alpha/pi
  double OneLoopVertexRho(double q2);
  //! rho^(2) = [F2/F1]^(2), in (alpha/pi)^2
  double TwoLoopVertexRho(double q2);

  //! T1/T0 and T2/T0 of the Pauli term (see the file comment)
  struct Pauli_Factors {
    double R1 = 0., R2 = 0.;
  };

  /*!
    R1, R2 for the line p -> pp of mass M, exchanging one photon with the
    point-like line k -> kp of mass m (Gordon decomposition, checked against
    explicit Dirac traces in yfs-helpers/two-loop-muon-vertex/pauli_factors.py).
  */
  Pauli_Factors PauliFactors(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &pp,
                             const ATOOLS::Vec4D &k, const ATOOLS::Vec4D &kp,
                             double M, double m);

  /*!
    X of the file comment, in (alpha/pi)^2, for the line p -> pp of mass M.
    The Pauli terms need the other line k -> kp (mass m); with withpauli false
    only [ln F1]^(2) is returned.
  */
  double CrossedLineTwoLoop(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &pp, double M,
                            bool withpauli, const ATOOLS::Vec4D &k,
                            const ATOOLS::Vec4D &kp, double m);

}

#endif
