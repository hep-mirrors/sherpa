/*!
  \file Acraic_LoopHelicity.C

  YFS: ACRAIC_Virtual: helicity. The loop provider's per-helicity virtual
  factor delta^h (YFS/NLO/NLO_LoopHelicity.C) on ACRAIC's coherent amplitude,

      rho_V = 1/4 sum_h | A_1^h + delta^h A_0^h |^2 ,

  which is ACRAIC_Virtual: external for delta^h = v/2. delta^h is computed
  once per event at the Born point and applied to every partition.
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "YFS/Main/Dyson_VP.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/MathTools.H"

#include <iomanip>

using namespace YFS;
using namespace ATOOLS;


/*
  Maps the provider's delta (physical labels hel, trees a0) to ACRAIC's d[f],
  unfilled entries keep vhalf. The leg-flip map is derived from the |a0|
  pattern against Comix's Born over the first ACRAIC: VIRT_HEL_NCAL events
  (ACRAIC: VIRT_HEL_MASK >= 0 overrides).
*/
bool Acraic_Base::LoopHelicityToAcraic(const std::vector<std::vector<int> > &hel,
                                   const std::vector<Complex> &a0,
                                   const std::vector<Complex> &delta,
                                   const Vec4D_Vector &p, Complex vhalf,
                                   std::vector<Complex> &d)
{
  const int nh(Amplitude::NHel());
  d.assign(nh, vhalf);
  if (hel.empty() || hel.size() != a0.size() || hel.size() != delta.size())
    return false;
  const size_t nl(m_flavs.size());
  // Provider labels packed as an ACRAIC index: bit Bit(i) set for helicity -1.
  // Legs with one spin state have no bit and must carry 0.
  std::vector<int> base(hel.size(), 0);
  for (size_t c(0); c < hel.size(); ++c) {
    if (hel[c].size() != nl) return false;
    int b(0);
    for (size_t i(0); i < nl; ++i) {
      const int bit(Amplitude::Bit(i));
      if (bit < 0) { if (hel[c][i] != 0) return false; continue; }
      if (hel[c][i] == 0) return false;
      if (hel[c][i] < 0) b = AddLeg(b, bit);
    }
    if (b >= nh) return false;
    base[c] = b;
  }
  if (m_comixflip < 0) return false;
  static const int maskset(Settings::GetMainSettings()["ACRAIC"]["VIRT_HEL_MASK"]
                           .SetDefault(-1).Get<int>());
  static const int ncal(Settings::GetMainSettings()["ACRAIC"]["VIRT_HEL_NCAL"]
                        .SetDefault(3).Get<int>());
  if (maskset >= 0) m_loopmask = maskset;
  else if (m_loopcal < ncal) {
    // At the Born only a leg relabelling and a normalisation can differ, so
    // match the unit-normalised moduli.
    Amplitude C;
    if (!ComixBornAmplitude(p, C)) { ++m_lh_mapbad; return m_loopmask >= 0 ? true : false; }
    double sa(0.), sc(0.);
    for (size_t c(0); c < a0.size(); ++c) sa += std::norm(a0[c]);
    for (int f = 0; f < nh; ++f) sc += std::norm(C.m_A[f]);
    if (!(sa > 0.) || !(sc > 0.)) { ++m_lh_mapbad; return false; }
    const double ra(1./sqrt(sa)), rc(1./sqrt(sc));
    int best(-1);
    double bmet(-1.), nmet(-1.);
    for (int m = 0; m < nh; ++m) {
      double met(0.);
      for (size_t c(0); c < a0.size(); ++c)
        met += sqr(std::abs(a0[c])*ra - std::abs(C.m_A[FlipLegs(base[c], m)])*rc);
      met = sqrt(met);
      if (best < 0 || met < bmet) { nmet = bmet; bmet = met; best = m; }
      else if (nmet < 0. || met < nmet) nmet = met;
    }
    // modulus ratio over live helicities: constant iff both use the same spin basis
    double rmin(0.), rmax(0.);
    int nlive(0);
    for (size_t c(0); c < a0.size(); ++c) {
      const Complex A(a0[c]), B(C.m_A[FlipLegs(base[c], best)]);
      if (std::abs(A) > 1e-3/ra && std::abs(B) > 1e-3/rc) {
        const double r(std::abs(A/B)*sqrt(sc/sa));
        if (nlive++ == 0) rmin = rmax = r;
        else { rmin = Min(rmin, r); rmax = Max(rmax, r); }
      }
    }
    const bool decisive(bmet < 1e-2 && (nmet < 0. || nmet > 10.*bmet));
    if (!decisive) {
      ++m_lh_mapbad;
      msg_Error()<<METHOD<<"(): provider -> Comix helicity map not decisive:"
                 <<" best mask "<<best<<" metric "<<bmet<<", runner-up "<<nmet
                 <<". Set ACRAIC: VIRT_HEL_MASK by hand.\n";
      if (m_loopmask < 0) return false;
    }
    else if (m_loopmask >= 0 && best != m_loopmask) {
      ++m_lh_mapbad;
      msg_Error()<<METHOD<<"(): provider -> Comix helicity map changed from "
                 <<m_loopmask<<" to "<<best<<"; keeping the first.\n";
    }
    else {
      if (m_loopmask < 0)
        msg_Info()<<"ACRAIC_Virtual: helicity - provider -> Comix leg-flip mask "
                  <<best<<" (metric "<<bmet<<", runner-up "<<nmet<<"),"
                  <<" provider -> ACRAIC index mask "
                  <<FlipLegs(best, ComixLegFlips(Amplitude::s_nlegs))
                  <<"; |A0_provider/A0_Comix| over "<<nlive
                  <<" live helicities in ["<<rmin<<", "<<rmax
                  <<"] (normalised to the sums).\n";
      m_loopmask = best;
      m_loopmet = bmet; m_loopnext = nmet;
      m_looprmin = rmin; m_looprmax = rmax;
    }
    ++m_loopcal;
  }
  if (m_loopmask < 0) return false;
  const int fm(ComixLegFlips(Amplitude::s_nlegs));
  for (size_t c(0); c < delta.size(); ++c)
    d[FlipLegs(FlipLegs(base[c], m_loopmask), fm)] = delta[c];
  return true;
}


int Acraic_Base::PhysicalHelicity(int f, size_t i) const
{
  const int bit(Amplitude::Bit(i));
  if (bit < 0) return 0;
  // ACRAIC index -> physical labels; without a map only f fbar is known
  int inv(0);
  if (m_loopmask >= 0 && m_comixflip >= 0)
    inv = FlipLegs(m_loopmask, ComixLegFlips(Amplitude::s_nlegs));
  else if (!m_ffbar) return 0;
  return LegHelicity(FlipLegs(f, inv), bit);
}


// 1/4 sum_f |A_1[f] + d[f] A_0[f]|^2 with A_0 = m_AmpExpo0 (partition-summed
// Born); d[f] = v/2 reproduces the external form.
double Acraic_Base::HelicityVirtualRho(const std::vector<Complex> &d,
                                     double P1, double P2) const
{
  const int nh(Amplitude::NHel());
  const bool pol(P1 != 0. || P2 != 0.);
  double sum(0.);
  for (int f = 0; f < nh && f < (int)d.size(); ++f) {
    double w(1.);
    if (pol) {
      const int l1(PhysicalHelicity(f, 0)), l2(PhysicalHelicity(f, 1));
      if (l1 == 0 || l2 == 0) return -1.;
      w = (1. + l1*P1)*(1. + l2*P2);
    }
    sum += w*std::norm(m_AmpExpo1.m_A[f] + d[f]*m_AmpExpo0.m_A[f]);
  }
  return sum/4.;
}


double Acraic_Base::A0InsertionRho(const std::vector<Complex> &e, double vhalf,
                                 bool vonbeta1) const
{
  const int nh(Amplitude::NHel());
  const Amplitude &ax(m_order == 2 ? m_AmpExpo2 : m_AmpExpo1);
  const Amplitude &av(vonbeta1 ? m_AmpExpo1 : m_AmpExpo0);
  double sum(0.);
  for (int f = 0; f < nh && f < (int)e.size(); ++f)
    sum += std::norm(ax.m_A[f] + vhalf*av.m_A[f] + e[f]*m_AmpExpo0.m_A[f]);
  return sum/4.;
}


std::vector<Complex> Acraic_Base::HelicityShares(const VP_Factors &vp) const
{
  const int nh(Amplitude::NHel());
  std::vector<Complex> r(nh, vp.rbar);
  if (Amplitude::s_nlegs != 4) return r;
  for (int f = 0; f < nh; ++f) {
    int l[4];
    bool known(true);
    for (size_t i(0); i < 4; ++i) {
      l[i] = PhysicalHelicity(f, i);
      if (l[i] == 0) known = false;
    }
    if (known) r[f] = vp.Share(l[0], l[1], l[2], l[3]);
  }
  return r;
}


std::vector<Complex> Acraic_Base::PhysicalLabelsToAcraic(const std::vector<std::vector<int> > &hel,
                                                     const std::vector<Complex> &val,
                                                     const Complex &fallback) const
{
  const int nh(Amplitude::NHel());
  const size_t nl(m_flavs.size());
  std::vector<Complex> r(nh, fallback);
  for (int f = 0; f < nh; ++f) {
    std::vector<int> l(nl, 0);
    bool known(true);
    for (size_t i(0); i < nl; ++i) {
      l[i] = PhysicalHelicity(f, i);
      if (l[i] == 0 && Amplitude::Bit(i) >= 0) known = false;
    }
    if (!known) continue;
    for (size_t c(0); c < hel.size() && c < val.size(); ++c)
      if (hel[c] == l) { r[f] = val[c]; break; }
  }
  return r;
}


double Acraic_Base::A0WeightedTwoRe(const std::vector<Complex> &e) const
{
  double num(0.), den(0.);
  for (int f = 0; f < Amplitude::NHel() && f < (int)e.size(); ++f) {
    num += 2.*e[f].real()*std::norm(m_AmpExpo0.m_A[f]);
    den += std::norm(m_AmpExpo0.m_A[f]);
  }
  return den > 0. ? num/den : 0.;
}


void Acraic_Base::HelicityVirtualSplit(const std::vector<Complex> &d, double vhalf,
                                     double split[3]) const
{
  split[0] = split[1] = split[2] = 0.;
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh && f < (int)d.size(); ++f) {
    const Complex D(d[f] - vhalf), a0(m_AmpExpo0.m_A[f]);
    const Complex b(m_AmpExpo1.m_A[f] + vhalf*a0);
    split[0] += 2.*std::real(std::conj(b)*D*a0);
    split[1] += sqr(D.real())*std::norm(a0);
    split[2] += sqr(D.imag())*std::norm(a0);
  }
  for (int i(0); i < 3; ++i) split[i] /= 4.;
}


void Acraic_Base::AccumulateLoopHelicity(bool ok, double addhel, double addext,
                                       double reg, const double split[3])
{
  m_lh_regmax = Max(m_lh_regmax, reg);
  if (!ok) { ++m_lh_fail; return; }
  ++m_lh_n;
  for (int i(0); i < 3; ++i) m_lh_split[i] += split[i];
  m_lh_sumhel += addhel;
  m_lh_sumext += addext;
  m_lh_sumd += addhel - addext;
  m_lh_sumd2 += sqr(addhel - addext);
  m_lh_dmax = Max(m_lh_dmax, std::abs(addhel - addext));
}


void Acraic_Base::ReportLoopHelicity()
{
  // same decision on every rank, so the reductions below stay matched
  if (!m_useacraic || m_acraicvirtsrc != acraicvirt::helicity) return;
  double buf[10] = {(double)m_lh_n, (double)m_lh_fail, (double)m_lh_mapbad,
                    m_lh_sumhel, m_lh_sumext, m_lh_sumd, m_lh_sumd2,
                    m_lh_split[0], m_lh_split[1], m_lh_split[2]};
  double mx[2] = {m_lh_regmax, m_lh_dmax};
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 10, MPI_DOUBLE, MPI_SUM);
    mpi->Allreduce(mx, 2, MPI_DOUBLE, MPI_MAX);
  }
#endif
  if (buf[0] + buf[1] <= 0.) return;
  const double n(buf[0] > 0. ? buf[0] : 1.);
  const double dm(buf[5]/n), drms(sqrt(Max(0., buf[6]/n - dm*dm)));
  msg_Out()<<std::setprecision(6)
           <<"ACRAIC rho_V (ACRAIC_Virtual: helicity): "<<(long)buf[0]
           <<" events helicity-resolved, "<<(long)buf[1]<<" at delta = v/2,"
           <<" "<<(long)buf[2]<<" map calibration problems; map mask "
           <<m_loopmask<<" (metric "<<m_loopmet<<", runner-up "<<m_loopnext<<")\n"
           <<"  <(rho_V - rho_1)/rho_crude>: helicity "<<buf[3]/n
           <<"  external "<<buf[4]/n<<"  difference mean "<<dm
           <<" rms "<<drms<<" max |.| "<<mx[1]<<"\n"
           <<"  split of the difference: linear in (delta - v/2) "<<buf[7]/n
           <<", (Re)^2 "<<buf[8]/n<<", (Im)^2 "<<buf[9]/n<<"\n"
           <<"  regression max |rho_V(delta=v/2) - external| / rho_1 = "
           <<mx[0]<<std::endl;
}
