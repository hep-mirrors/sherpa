/*!
  \file Beta0_2_VPInsertion.C

  The dispersive lepton-loop insertion in a space-like vertex; see Beta0_2_VPInsertion.H.
*/

#include "YFS/Main/Beta0_2_VPInsertion.H"
#include "YFS/Main/Beta0_2_Vertex.H"
#include "YFS/Main/Beta0_2_VP_Table.H"
#include "YFS/Main/Beta0_2_VPBox_Table.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"

#include <cmath>

using namespace YFS;
using namespace ATOOLS;

namespace {

  //! Catmull-Rom interpolation of a column on the uniform grid l0 + i dl (inside it)
  double Interpolate(const double *v, int n, double l0, double dl, double l)
  {
    const double u(Max(0., Min(double(n - 1), (l - l0)/dl)));
    const int i(Min(n - 2, int(u)));
    const double t(u - i);
    const double p0(v[Max(0, i - 1)]), p1(v[i]), p2(v[i + 1]), p3(v[Min(n - 1, i + 2)]);
    return p1 + 0.5*t*(p2 - p0 + t*(2.*p0 - 5.*p1 + 4.*p2 - p3
                                    + t*(3.*(p1 - p2) + p3 - p0)));
  }

  bool SameMass(double a, double b) { return std::abs(a - b) <= 1e-6*b; }

}


bool YFS::VPInsertionMassesMatch()
{
  using namespace B02VPTable;
  static const bool ok(SameMass(Flavour(kf_e).Mass(true), s_M_e)
                       && SameMass(Flavour(kf_mu).Mass(true), s_M_mu)
                       && SameMass(Flavour(kf_tau).Mass(true), s_M_tau));
  static bool told(false);
  if (!ok && !told) {
    told = true;
    msg_Error()<<"YFS: BETA0_2_VP_INSERTIONS: the table is built for m_e = "<<s_M_e
               <<", m_mu = "<<s_M_mu<<", m_tau = "<<s_M_tau<<" GeV (PDG), the run has "
               <<Flavour(kf_e).Mass(true)<<", "<<Flavour(kf_mu).Mass(true)<<", "
               <<Flavour(kf_tau).Mass(true)<<": no VP insertion.\n";
  }
  return ok;
}


double YFS::VPInsertionVertex(const Vec4D &p, const Vec4D &pp, double M,
                              const Vec4D &k, const Vec4D &kp, double m)
{
  using namespace B02VPTable;
  if (!VPInsertionMassesMatch()) return 0.;
  const double q2(-(p - pp).Abs2());
  if (!(q2 > 0.)) return 0.;
  const bool eline(SameMass(M, s_M_e)), muline(SameMass(M, s_M_mu));
  if (!eline && !muline) return 0.;
  const int n(eline ? s_n_e : s_n_mu);
  const double l0(eline ? s_l0_e : s_l0_mu), dl(eline ? s_dl_e : s_dl_mu);
  const double l(log(q2/(M*M)));
  if (l < l0 || l > l0 + (n - 1)*dl) {
    static bool told(false);
    if (!told) {
      told = true;
      msg_Error()<<"YFS: BETA0_2_VP_INSERTIONS: ln(q2/M^2) = "<<l<<" outside the table ["
                 <<l0<<", "<<l0 + (n - 1)*dl<<"]: no VP insertion there (reported once).\n";
    }
    return 0.;
  }
  double f1(0.), f2(0.);
  for (int loop(0); loop < 3; ++loop) {
    f1 += Interpolate(eline ? s_F1_e[loop] : s_F1_mu[loop], n, l0, dl, l);
    f2 += Interpolate(eline ? s_F2_e[loop] : s_F2_mu[loop], n, l0, dl, l);
  }
  return f1 + f2*PauliFactors(p, pp, k, kp, M, m).R1;
}


double YFS::VPInsertionBox(const Flavour_Vector &fl, double s, double t)
{
  using namespace B02VPBoxTable;
  if (!VPInsertionMassesMatch() || fl.size() != 4) return 0.;
  const kf_code a(fl[0].Kfcode()), b(fl[1].Kfcode());
  if (!((a == kf_mu && b == kf_e) || (a == kf_e && b == kf_mu))) return 0.;
  if (!(t < 0.)) return 0.;
  int ib(0);
  for (int i(1); i < s_nbeam; ++i)
    if (std::abs(s - s_s[i]) < std::abs(s - s_s[ib])) ib = i;
  // the grid ends just below |t|max of its beam energy; allow rounding above it, nothing more
  const double l(log(-t)), lmax(s_l0[ib] + (s_n - 1)*s_dl[ib]);
  if (l < s_l0[ib] || l > lmax + 1e-6) {
    static bool told(false);
    if (!told) {
      told = true;
      msg_Error()<<"YFS: BETA0_2_VP_INSERTIONS vertex_box: ln|t| = "<<l<<" outside the box table ["
                 <<s_l0[ib]<<", "<<lmax<<"] at s = "<<s<<": no box insertion there (reported once).\n";
    }
    return 0.;
  }
  return fl[0].Charge()*fl[1].Charge()*Interpolate(s_N[ib], s_n, s_l0[ib], s_dl[ib], Min(l, lmax));
}
