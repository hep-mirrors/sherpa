/*!
  \file Acraic_RealVirtualHelicity.C

  ACRAIC: REAL_VIRTUAL helicity. The non-factorisable O(alpha^2)
  real-virtual resolved by helicity: the loop provider's per-helicity
  one-loop amplitudes of the (n+1)-point process (Recola, read after the
  RV_MODE remainder loop call, NLO_Base::BuildRealVirtualHelicityFactors)
  put on photon j's M_1 terms of ACRAIC's coherent amplitude,

      rho += 1/4 sum_h sum_j 2 Re[ conj(A^h) r_j^h M_1j^h ] ,

      r_j^h = delta_{n+1,j}^h - delta_B^h ,

  A = A_1 (A_2 at ORDER 2), so that M_1j carries the (n+1)-point's own loop
  over tree per helicity and its eikonal x Born subtraction the Born's. This
  is the helicity-resolved form of REAL_VIRTUAL averaged, which is the case
  r_j^h = dv_j/2 on every h (RealVirtualRemainderRho); like it, linear in
  the remainder (O(alpha^2)), with the factorisable v/2 placed elsewhere.

  In the one-loop N-photon language (NOTES-acraic-oneloop-nphoton.md) this
  is the n = 1 member: beta_1^[1](k) = [delta_1 - delta_0] M_1 + delta_0
  beta_1^[0], the second term being REAL_VIRTUAL factorisable.

  Partitions: photon j's M_1 is summed over the partitions that carry it
  (m_realphotM1[j], per-partition propagator shifts), while r_j^h comes from
  ONE (n+1)-point per photon, as the Born's delta^h comes from one Born
  point (the "ratio" placement of the notes, sec. 5.2).
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "YFS/NLO/NLO_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/MathTools.H"

#include <iomanip>

using namespace YFS;
using namespace ATOOLS;

namespace {

  /*!
    The provider's (n+1)-leg labels as a packed index: the fermion legs as
    LoopHelicityToAcraic packs them (bit Bit(i) set for helicity -1), the
    photon in bit nlg (set for -1). -1 if a label has no place.
  */
  int PackedLabels(const std::vector<int> &h, size_t nl, int nlg)
  {
    if (h.size() != nl + 1) return -1;
    int b(0);
    for (size_t i(0); i < nl; ++i) {
      const int bit(Amplitude::Bit(i));
      if (bit < 0) { if (h[i] != 0) return -1; continue; }
      if (h[i] == 0) return -1;
      if (h[i] < 0) b |= (1 << bit);
    }
    if (h[nl] == 0) return -1;
    if (h[nl] < 0) b |= (1 << nlg);
    return b;
  }

}


bool Acraic_Base::RealVirtualHelicityToAcraic(const RVHelicityFactors &rv,
                                              int phohel, Complex dvhalf,
                                              std::vector<Complex> &r)
{
  const int nh(Amplitude::NHel()), nlg(Amplitude::s_nlegs);
  r.assign(nh, dvhalf);
  if (!rv.ok || rv.hel.empty() || rv.hel.size() != rv.rem.size()
      || rv.hel.size() != rv.a0.size() || phohel == 0) return false;
  const size_t nl(m_flavs.size());
  std::vector<int> base(rv.hel.size(), 0);
  for (size_t c(0); c < rv.hel.size(); ++c) {
    base[c] = PackedLabels(rv.hel[c], nl, nlg);
    if (base[c] < 0 || base[c] >= 2*nh) return false;
  }
  if (m_comixflip < 0) return false;
  static const int maskset(Settings::GetMainSettings()["ACRAIC"]["RV_HEL_MASK"]
                           .SetDefault(-1).Get<int>());
  static const int ncal(Settings::GetMainSettings()["ACRAIC"]["RV_HEL_NCAL"]
                        .SetDefault(3).Get<int>());
  if (maskset >= 0) m_rvmask = maskset;
  else if (m_rvcal < ncal) {
    /*
      Calibration at the (n+1)-point: Comix's one-photon amplitude at the
      same lab point, both photon helicities, already in ACRAIC's index
      (ComixRealAt applies m_comixflip), against the provider's tree over
      every relabelling of the n+1 legs. Both patterns normalised to unit
      sum. The moduli of definite-helicity amplitudes carry no phase
      convention, for the photon included.
    */
    Amplitude Cp, Cm;
    if (rv.plab.size() != nl + 1
        || !ComixRealAt(rv.plab, +1, Cp, -1.) || !ComixRealAt(rv.plab, -1, Cm, -1.)) {
      ++m_rvh_mapbad;
      return m_rvmask >= 0;
    }
    double sa(0.), sc(0.);
    for (size_t c(0); c < rv.a0.size(); ++c) sa += std::norm(rv.a0[c]);
    for (int f = 0; f < nh; ++f) sc += std::norm(Cp.m_A[f]) + std::norm(Cm.m_A[f]);
    if (!(sa > 0.) || !(sc > 0.)) { ++m_rvh_mapbad; return false; }
    const double ra(1./sqrt(sa)), rc(1./sqrt(sc));
    auto comix = [&](int idx) -> double {
      const Amplitude &C((idx >> nlg) & 1 ? Cm : Cp);
      return std::abs(C.m_A[idx & (nh - 1)]);
    };
    int best(-1);
    double bmet(-1.), nmet(-1.);
    for (int m = 0; m < 2*nh; ++m) {
      double met(0.);
      for (size_t c(0); c < rv.a0.size(); ++c)
        met += sqr(std::abs(rv.a0[c])*ra - comix(base[c] ^ m)*rc);
      met = sqrt(met);
      if (best < 0 || met < bmet) { nmet = bmet; bmet = met; best = m; }
      else if (nmet < 0. || met < nmet) nmet = met;
    }
    const bool decisive(bmet < 1e-2 && (nmet < 0. || nmet > 10.*bmet));
    const int bornmap(m_loopmask >= 0 ? (m_loopmask ^ (m_comixflip & (nh - 1))) : -1);
    if (!decisive) {
      ++m_rvh_mapbad;
      msg_Error()<<METHOD<<"(): provider -> ACRAIC (n+1)-point helicity map not"
                 <<" decisive: best mask "<<best<<" metric "<<bmet
                 <<", runner-up "<<nmet<<". Set ACRAIC: RV_HEL_MASK by hand.\n";
      if (m_rvmask < 0) return false;
    }
    else if (m_rvmask >= 0 && best != m_rvmask) {
      ++m_rvh_mapbad;
      msg_Error()<<METHOD<<"(): provider -> ACRAIC (n+1)-point helicity map"
                 <<" changed from "<<m_rvmask<<" to "<<best<<"; keeping the first.\n";
    }
    else {
      if (m_rvmask < 0) {
        msg_Info()<<"ACRAIC: REAL_VIRTUAL helicity - provider -> ACRAIC (n+1)-point"
                  <<" map mask "<<best<<" (fermion bits "<<(best & (nh - 1))
                  <<", photon flip "<<((best >> nlg) & 1)<<"; metric "<<bmet
                  <<", runner-up "<<nmet<<"); Born map "<<bornmap<<".\n";
        if (bornmap >= 0 && (best & (nh - 1)) != bornmap) {
          ++m_rvh_mapbad;
          msg_Error()<<METHOD<<"(): the fermion part of the (n+1)-point map, "
                     <<(best & (nh - 1))<<", differs from the Born map "
                     <<bornmap<<".\n";
        }
      }
      m_rvmask = best;
      m_rvmet = bmet; m_rvnext = nmet;
    }
    ++m_rvcal;
  }
  if (m_rvmask < 0) return false;
  const int want(phohel > 0 ? 0 : 1);
  size_t nset(0);
  for (size_t c(0); c < rv.rem.size(); ++c) {
    const int idx(base[c] ^ m_rvmask);
    if (((idx >> nlg) & 1) != want) continue;
    r[idx & (nh - 1)] = rv.rem[c];
    ++nset;
  }
  return nset > 0;
}


double Acraic_Base::HelicityRealVirtualRemainderRho
  (const std::vector<std::vector<Complex> > &r) const
{
  const Amplitude &A(m_order == 2 ? m_AmpExpo2 : m_AmpExpo1);
  const int nh(Amplitude::NHel());
  double sum(0.);
  for (size_t j(0); j < r.size() && j < m_realphotM1.size(); ++j)
    for (int f = 0; f < nh && f < (int)r[j].size(); ++f)
      sum += 2.*std::real(conj(A.m_A[f]) * r[j][f] * m_realphotM1[j].m_A[f]);
  return sum/4.;
}


void Acraic_Base::AccumulateRealVirtualHelicity(bool anyhel, double addhel,
                                                double addavg)
{
  if (!anyhel) { ++m_rvh_avg; return; }
  ++m_rvh_n;
  m_rvh_sumhel += addhel;
  m_rvh_sumavg += addavg;
  m_rvh_sumd += addhel - addavg;
  m_rvh_sumd2 += sqr(addhel - addavg);
  m_rvh_dmax = Max(m_rvh_dmax, std::abs(addhel - addavg));
}


void Acraic_Base::ReportRealVirtualHelicity()
{
  // same decision on every rank, so the reductions below stay matched
  if (RealVirtualMode() != acraicrv::helicity) return;
  double buf[8] = {(double)m_rvh_n, (double)m_rvh_avg, (double)m_rvh_mapbad,
                   m_rvh_sumhel, m_rvh_sumavg, m_rvh_sumd, m_rvh_sumd2, 0.};
  double mx[1] = {m_rvh_dmax};
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 8, MPI_DOUBLE, MPI_SUM);
    mpi->Allreduce(mx, 1, MPI_DOUBLE, MPI_MAX);
  }
#endif
  if (buf[0] + buf[1] <= 0.) return;
  const double n(buf[0] > 0. ? buf[0] : 1.);
  const double dm(buf[5]/n), drms(sqrt(Max(0., buf[6]/n - dm*dm)));
  msg_Out()<<std::setprecision(6)
           <<"ACRAIC: REAL_VIRTUAL helicity: "<<(long)buf[0]
           <<" events with a helicity-resolved remainder, "<<(long)buf[1]
           <<" with none (averaged or no RV photon), "<<(long)buf[2]
           <<" map calibration problems; (n+1)-point map mask "<<m_rvmask
           <<" (metric "<<m_rvmet<<", runner-up "<<m_rvnext<<")\n"
           <<"  <added/rho_crude> on those events: helicity RV "
           <<buf[3]/n<<"  averaged RV "<<buf[4]/n<<"  difference mean "<<dm
           <<" rms "<<drms<<" max |.| "<<mx[0]<<std::endl;
}
