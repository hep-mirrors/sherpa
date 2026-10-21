/*
  The YFS_Handler's diagnostics and probes: per-event print-outs switched on
  by YFS / ACRAIC settings (WEIGHT_PROBE, LO_PROBE, RR_PROBE, ...). None of
  them changes a weight.
*/
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <sstream>
#include <iomanip>
#include "YFS/Main/YFS_Handler.H"
#include "YFS/Main/ISR.H"

using namespace std;
using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;
using namespace PHASIC;
using namespace METOOLS;

namespace {
  /*
    The eikonal of a final-state pair without couplings,
    -(p1/(p1.k) - p2/(p2.k))^2 = 2 p1.p2/(p1.k p2.k) - m1^2/(p1.k)^2 - m2^2/(p2.k)^2.
  */
  double PairEikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2)
  {
    return -(p1/(p1*k) - p2/(p2*k)).Abs2();
  }

  /*
    YFS: FSR_DENSITY_GATE (diagnostic, default 0 = off; N > 0 prints the
    first N events with FSR photons). Which eikonal is the density the FSR
    generator actually produced its photons with?

    FSR::GenerateAngles samples the angle from the crude f-bar, which uses the
    PRE-emission velocities; FSR::F then multiplies every photon by the mass
    weight m_f/f-bar, and m_f (KKMC's KarFin dist1) is built from the
    POST-emission pair in its own rest frame. The weight is in the event
    weight (YFS: WEIGHT_MODE Full, the default), so the weighted density is
    m_f: 4 m_f/w'^2 = S(post legs; k), w' the photon energy in the
    post-emission pair frame. One line per FSR photon:
      gpost = 4 m_f/(w'^2 S(post))  - the generator against the eikonal on
              the physical legs (ACRAIC's crude, up to the coupling),
      gpre  = 4 m_f/(w'^2 S(pre))   - against the eikonal on the
              pre-emission legs (YFS.NLO's FSR crude, REAL_SUB_EIK 8),
      pfyy  = (q'+K)^2/q'^2 * FSR::m_yy - ACRAIC's pseudo-flux against the
              generator's rescaling (1 if the flux is the generator's),
      mw    = FSR's mass weight / (m_f/f-bar) (1 if the weight is m_f/f-bar),
      resid = |Q_pre - q1' - q2' - K| (the pre-emission legs, the physical
              legs and the photons in one frame),
    with thk = the photon angle to the nearer physical lepton in units of
    m/E of that lepton and x = 2 k0/sqrt(s).
  */
  void PrintFSRDensityGate(const FSR &fsr, YFS::Dipole &D,
                           const Vec4D_Vector &prefsr, double sqrts,
                           size_t nisr, double m1, double m2)
  {
    const int l(D.Left()), r(D.Right());
    if (l < 0 || r < 0 || l >= (int)prefsr.size() || r >= (int)prefsr.size()) return;
    const Vec4D q1(D.GetNewMomenta(0)), q2(D.GetNewMomenta(1));
    const Vec4D p1(prefsr[l]), p2(prefsr[r]);
    const Vec4D_Vector k(D.GetPhotons());
    Vec4D K;
    for (const Vec4D &g : k) K += g;
    const Vec4D Qp(q1 + q2), res(p1 + p2 - Qp - K);
    const double sq(Qp.Abs2());
    const double pfyy(sq > 0. ? (Qp + K).Abs2()/sq*fsr.m_yy : -1.);
    const double resid(sqrt(sqr(res[0]) + sqr(res[1]) + sqr(res[2]) + sqr(res[3])));
    // the generation entries, by their energy in the post-emission frame
    std::vector<double> wgen;
    for (size_t i(0); i < fsr.m_k0.size(); ++i)
      wgen.push_back(fsr.m_k0[i]*fsr.m_xfact*fsr.m_scalek);
    for (const Vec4D &g : k) {
      const double w((g*Qp)/sqrt(sq));
      size_t best(wgen.size()); double dmin(1e-6);
      for (size_t i(0); i < wgen.size(); ++i)
        if (std::abs(wgen[i]/w - 1.) < dmin) { dmin = std::abs(wgen[i]/w - 1.); best = i; }
      if (best == wgen.size() || best >= fsr.m_dist1.size()) {
        std::cerr<<"@@@ FSRGATE unmatched w="<<w<<std::endl;
        continue;
      }
      const double mf(fsr.m_dist1[best]), fb(fsr.m_dist2[best]);
      const double spost(PairEikonal(g, q1, q2)), spre(PairEikonal(g, p1, p2));
      const double gpost(4.*mf/(w*w*spost)), gpre(4.*mf/(w*w*spre));
      const double mw(best < fsr.m_MassWls.size() && fb > 0. && mf > 0.
                      ? fsr.m_MassWls[best]/(mf/fb) : -1.);
      const double c1(Vec3D(g)*Vec3D(q1)/(Vec3D(g).Abs()*Vec3D(q1).Abs()));
      const double c2(Vec3D(g)*Vec3D(q2)/(Vec3D(g).Abs()*Vec3D(q2).Abs()));
      const bool one(c1 > c2);
      const Vec4D &qn(one ? q1 : q2);
      const double th(acos(Min(1., one ? c1 : c2)));
      const double thk(th/((one ? m1 : m2)/qn[0]));
      std::cerr<<std::setprecision(8)<<"@@@ FSRGATE nisr="<<nisr<<" nfsr="<<k.size()
               <<" x="<<2.*g[0]/sqrts<<" thk="<<thk<<" gpost="<<gpost
               <<" gpre="<<gpre<<" pfyy="<<pfyy<<" mw="<<mw<<" resid="<<resid
               <<" spost="<<spost<<std::endl;
    }
  }
}


void YFS_Handler::ProbeFSRDensity(YFS::Dipole &D)
{
  static const int gate(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["FSR_DENSITY_GATE"].SetDefault(0).Get<int>());
  static int ngate(0);
  if (gate > 0 && ngate < gate && !D.GetPhotons().empty()) {
    ++ngate;
    PrintFSRDensityGate(*p_fsr, D, m_ev.m_reallab, sqrt(m_s),
                        m_ev.m_ISRPhotons.size(), m_flavs[D.Left()].Mass(),
                        m_flavs[D.Right()].Mass());
  }
}
