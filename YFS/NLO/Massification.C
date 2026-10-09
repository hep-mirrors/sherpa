#include "YFS/NLO/Massification.H"
#include "YFS/Main/YFS_Options.H"

#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <cmath>

using namespace YFS;
using namespace ATOOLS;

namespace {
  // setup bound on 4 m^2/s of a massified species: the e+e- -> mu mu residual
  const double s_max_4m2_over_s = 1e-3;
  // |t| (or |s_ij|) below this many m^2: the forward region where the
  const double s_forward_m2 = 100.;
  // the double pole is analytic in both: a relative mismatch above this is a
  const double s_double_pole_tolerance = 1e-8;
}

bool Massification::Enabled()
{
  return Mode() != virtmassif::off;
}

virtmassif::code Massification::Mode()
{
  static const virtmassif::code mode(Settings::GetMainSettings()["YFS"]
      ["VIRTUAL_MASSIFICATION"].SetDefault(virtmassif::off).Get<virtmassif::code>());
  return mode;
}

namespace {
  //! whether YFS: VIRTUAL_MASSIFICATION makes this species massless
  bool Massified(const kf_code kf, const virtmassif::code mode)
  {
    switch (mode) {
    case virtmassif::off:            return false;
    case virtmassif::electrons:      return kf == kf_e;
    case virtmassif::leptons:        return kf == kf_e || kf == kf_mu;
    case virtmassif::light_fermions: return kf == kf_e || kf == kf_mu ||
                                            kf == kf_u || kf == kf_d || kf == kf_s;
    }
    return false;
  }
}

Massification::Massification(const Flavour_Vector &flavs, const std::string &what) :
  m_what(what), m_flavs(flavs), m_c(0.), m_n(0), m_nforward(0), m_ndeadcone(0)
{
  const double s(sqr(rpa->gen.Ecms()));
  for (size_t i(0); i < flavs.size(); ++i) {
    const kf_code kf(flavs[i].Kfcode());
    m_provmass.push_back(flavs[i].Mass());
    if (flavs[i].IsPhoton()) m_photons.push_back(i);
    if (!Massified(kf, Mode())) continue;
    const double m(flavs[i].Mass());
    if (!(m > 0.))
      THROW(fatal_error, "YFS: VIRTUAL_MASSIFICATION needs Massive: true for "
                         + flavs[i].IDName() + ": the massless provider value"
                         " is rebuilt into the massive one");
    if (4.*m*m/s > s_max_4m2_over_s)
      THROW(fatal_error, "YFS: VIRTUAL_MASSIFICATION: 4 m^2/s = "
                         + ToString(4.*m*m/s) + " for " + flavs[i].IDName()
                         + " exceeds " + ToString(s_max_4m2_over_s)
                         + " - too close to the pair threshold (Coulomb)");
    m_legs.push_back(i);
    m_provmass[i] = 0.;
    PHASIC::Virtual_ME2_Base::AddMasslessLoopSpecies(kf);
  }
}

void Massification::Report()
{
  if (m_legs.empty()) return;
  long buf[4] = {m_n, m_nforward, m_ndeadcone, m_ne2bad};
  double worst(m_e2worst);
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 4, MPI_LONG, MPI_SUM);
    mpi->Allreduce(&worst, 1, MPI_DOUBLE, MPI_MAX);
  }
#endif
  msg_Out()<<"YFS: VIRTUAL_MASSIFICATION ("<<m_what<<"): "<<buf[0]
           <<" evaluations; forward (|s_ij| < "<<s_forward_m2<<" m^2): "
           <<buf[1]<<"; photon inside m/E of a massless lepton: "<<buf[2]
           <<"; 1/eps^2 not cancelled: "<<buf[3]
           <<(buf[3] ? " (worst relative "+ToString(worst)+")" : "")
           <<std::endl;
  if (buf[3]) msg_Error()<<"YFS: VIRTUAL_MASSIFICATION ("<<m_what<<"): the "
                         <<"double pole did not cancel on "<<buf[3]
                         <<" evaluations; the massified virtual is not trustworthy"
                         <<std::endl;
}

void Massification::Attach(const PHASIC::Virtual_ME2_Base *me)
{
  if (m_legs.empty()) return;
  if (!me || !me->MasslessLegConstant(m_c))
    THROW(fatal_error, "YFS: VIRTUAL_MASSIFICATION: the " + m_what +
          " loop provider keeps the lepton masses (only OpenLoops implements"
          " it; Recola treats fermions below RECOLA_LIGHT_FERMION_THRESHOLD"
          " as light, mass in the logs only - the same at one loop)");
  msg_Debugging()<<"YFS: VIRTUAL_MASSIFICATION ("<<m_what<<"): "<<m_legs.size()
            <<" massless legs, provider constant "<<m_c<<" per leg\n";
}

Vec4D_Vector Massification::MasslessImage(const Vec4D_Vector &p) const
{
  if (m_legs.empty()) return p;
  const Vec4D P(p[0] + p[1]);
  const double rs(sqrt(P.Abs2()));
  Poincare cms(P);
  Vec4D_Vector q(p);
  for (Vec4D &v : q) cms.Boost(v);
  // incoming: two-body kinematics with the provider masses along the
  // original direction of the first beam
  const double m0(m_provmass[0]), m1(m_provmass[1]);
  const double E0((rs*rs + m0*m0 - m1*m1)/(2.*rs));
  const double pz(sqrt(Max(0., E0*E0 - m0*m0)));
  const Vec3D n(Vec3D(q[0])/Vec3D(q[0]).Abs());
  q[0] = Vec4D(E0, pz*n);
  q[1] = Vec4D(rs - E0, -pz*n);
  // outgoing: one common rescale xi of the three-momenta, by Newton from
  // xi = 1 (the energy sum is convex and increasing in xi)
  double xi(1.);
  for (int it(0); it < 100; ++it) {
    double f(-rs), df(0.);
    for (size_t i(2); i < q.size(); ++i) {
      const double a(Vec3D(q[i]).Sqr()), e(sqrt(xi*xi*a + sqr(m_provmass[i])));
      f += e;
      if (e > 0.) df += xi*a/e;
    }
    if (!(df > 0.)) break;
    const double step(f/df);
    xi -= step;
    if (std::abs(step) < 1e-15*xi) break;
  }
  for (size_t i(2); i < q.size(); ++i) {
    const Vec3D v(xi*Vec3D(q[i]));
    q[i] = Vec4D(sqrt(v.Sqr() + sqr(m_provmass[i])), v);
  }
  for (Vec4D &v : q) cms.BoostBack(v);
  return q;
}

double Massification::Z(const double muir2) const
{
  // per leg Q^2 [L^2 + L + 4 + zeta2 + c] in alpha/4pi; /2 for alpha/2pi
  double z(0.);
  for (const size_t i : m_legs) {
    const double L(log(muir2/sqr(m_flavs[i].Mass())));
    z += sqr(m_flavs[i].Charge())*(L*L + L + 4. + M_PI*M_PI/6. + m_c);
  }
  return 0.5*z;
}

double Massification::ZPole1(const double muir2) const
{
  double z(0.);
  for (const size_t i : m_legs)
    z += sqr(m_flavs[i].Charge())*(2.*log(muir2/sqr(m_flavs[i].Mass())) + 1.);
  return 0.5*z;
}

double Massification::ZPole2() const
{
  double z(0.);
  for (const size_t i : m_legs) z += sqr(m_flavs[i].Charge());
  return z;
}

void Massification::CheckDoublePole(const double providere2)
{
  if (m_legs.empty()) return;
  const double z(ZPole2());
  const double rel(std::abs(providere2 + z)/Max(std::abs(providere2), std::abs(z)));
  if (!(rel < s_double_pole_tolerance)) {
    ++m_ne2bad;
    m_e2worst = Max(m_e2worst, IsBad(rel) ? 1. : rel);
    if (m_ne2bad <= 5)
      msg_Error()<<"YFS: VIRTUAL_MASSIFICATION ("<<m_what<<"): the provider's "
                 <<"1/eps^2 "<<providere2<<" is not cancelled by Z's "<<z
                 <<" (relative "<<rel<<"): its massless double pole or its "
                 <<"epsilon convention is not the one Z assumes"<<std::endl;
  }
}

void Massification::Count(const Vec4D_Vector &p)
{
  if (m_legs.empty()) return;
  ++m_n;
  bool forward(false), deadcone(false);
  for (size_t a(0); a < m_legs.size(); ++a)
    for (size_t b(a + 1); b < m_legs.size(); ++b) {
      const size_t i(m_legs[a]), j(m_legs[b]);
      // crossing: an incoming and an outgoing leg form a t-like invariant
      const bool cross((i < 2) != (j < 2));
      const double inv(cross ? (p[i] - p[j]).Abs2() : (p[i] + p[j]).Abs2());
      const double m2(Max(sqr(m_flavs[i].Mass()), sqr(m_flavs[j].Mass())));
      if (std::abs(inv) < s_forward_m2*m2) forward = true;
    }
  for (const size_t k : m_photons)
    for (const size_t i : m_legs) {
      const double th(p[i].Theta(p[k]));
      if (th < m_flavs[i].Mass()/p[i][0]) deadcone = true;
    }
  if (forward) ++m_nforward;
  if (deadcone) ++m_ndeadcone;
}
