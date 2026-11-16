/*!
  \file Beta0_2_Vertex.C

  The exact photonic two-loop vertex of a space-like line; see Beta0_2_Vertex.H.
*/

#include "YFS/Main/Beta0_2_Vertex.H"
#include "YFS/Main/Beta0_2.H"
#include "YFS/Main/Beta0_2_Vertex_Table.H"
#include "ATOOLS/Math/MathTools.H"

#include <cmath>

using namespace YFS;
using namespace ATOOLS;

namespace {

  //! Catmull-Rom interpolation of one table column at l = ln q2 (inside the grid)
  double Interpolate(const double *v, double l)
  {
    using namespace B02VertexTable;
    const double u(Max(0., Min(double(s_n - 1), (l - s_l0)/s_dl)));
    const int i(Min(s_n - 2, int(u)));
    const double t(u - i);
    const double p0(v[Max(0, i - 1)]), p1(v[i]), p2(v[i + 1]),
                 p3(v[Min(s_n - 1, i + 2)]);
    return p1 + 0.5*t*(p2 - p0 + t*(2.*p0 - 5.*p1 + 4.*p2 - p3
                                    + t*(3.*(p1 - p2) + p3 - p0)));
  }

  double TableLmax()
  {
    using namespace B02VertexTable;
    return s_l0 + (s_n - 1)*s_dl;
  }

  //! C L + K of hep-ph/0508127, the q2 -> infinity limit of [ln F1]^(2)
  double HighEnergyLnF1(double l)
  {
    return Beta02NLLCoefficient()*(l - 1.) + Beta02TwoLoopConstant();
  }

}


double YFS::TwoLoopVertexLnF1(double q2)
{
  using namespace B02VertexTable;
  if (!(q2 > 0.)) return 0.;
  const double l(log(q2));
  // linear in q2 at small q2 (F1 - 1 vanishes like q2)
  if (l < s_l0) return s_lnF1[0]*q2/exp(s_l0);
  const double lmax(TableLmax());
  if (l > lmax) return HighEnergyLnF1(l) + (s_lnF1[s_n - 1] - HighEnergyLnF1(lmax));
  return Interpolate(s_lnF1, l);
}


double YFS::OneLoopVertexRho(double q2)
{
  using namespace B02VertexTable;
  if (!(q2 > 0.)) return s_rho1[0];
  return Interpolate(s_rho1, log(q2));
}


double YFS::TwoLoopVertexRho(double q2)
{
  using namespace B02VertexTable;
  if (!(q2 > 0.)) return s_rho2[0];
  return Interpolate(s_rho2, log(q2));
}


Pauli_Factors YFS::PauliFactors(const Vec4D &p, const Vec4D &pp,
                                const Vec4D &k, const Vec4D &kp,
                                double M, double m)
{
  Pauli_Factors out;
  const double M2(M*M), m2(m*m);
  const double pk(p*k), ppkp(pp*kp), ppk(pp*k), pkp(p*kp), ppp(p*pp), kkp(k*kp);
  // T0/16: the contraction of the two point-like tensors
  const double t0(2.*(pk*ppkp + ppk*pkp) - 2.*kkp*(ppp - M2) - 2.*ppp*(kkp - m2)
                  + 4.*(kkp - m2)*(ppp - M2));
  if (!(dabs(t0) > 0.)) return out;
  // L_ab P^a P^b / 16 with P = p + pp
  const Vec4D P(p + pp);
  const double lpp(0.25*(2.*(k*P)*(kp*P) - P.Abs2()*(kkp - m2)));
  const double ratio(lpp/t0);
  out.R1 = 1. - 2.*ratio;
  out.R2 = 1. + ratio*((ppp + M2)/M2 - 4.);
  return out;
}


double YFS::CrossedLineTwoLoop(const Vec4D &p, const Vec4D &pp, double M,
                               bool withpauli, const Vec4D &k, const Vec4D &kp,
                               double m)
{
  const double q2(-(p - pp).Abs2()/(M*M));
  double x(TwoLoopVertexLnF1(q2));
  if (!withpauli) return x;
  const Pauli_Factors r(PauliFactors(p, pp, k, kp, M, m));
  const double rho1(OneLoopVertexRho(q2));
  x += TwoLoopVertexRho(q2)*r.R1 + 0.5*rho1*rho1*(r.R2 - 2.*r.R1*r.R1);
  return x;
}
