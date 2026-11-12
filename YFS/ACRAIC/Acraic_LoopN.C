/*!
  \file Acraic_LoopN.C

  ACRAIC: LOOP_BETA_MAX 2 - the one-loop two-photon residual. For a photon pair S = {j, l} the
  one-loop subset recursion
    beta_S^[1] = delta_S M_S - s_l beta_j^[1] - s_j beta_l^[1] - s_j s_l beta_{}^[1],
    beta_j^[1] = delta_j M_1j - s_j delta_0 B,   beta_{}^[1] = delta_0 B,
  minus what the factorisable virtual (delta_0 on every beta) and REAL_VIRTUAL
  helicity (rem_j = delta_j - delta_0 on photon j's M_1) already put in, is
    X_S = (delta_S - delta_0) M_S - s_l rem_j M_1j - s_j rem_l M_1l .
  Soft l: M_S -> s_l M_1j, delta_S -> delta_j, so X_S -> -s_j rem_l M_1l -> 0.
  The deltas come from one (n+2)-point per pair ("ratio" placement), so the
  partition sum factorises into P_S, Qa_S, Qb_S (filled in Acraic_Beta2.C).
  Linear in the remainders, O(alpha^3).
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "YFS/NLO/NLO_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <sstream>

#include <iomanip>

using namespace YFS;
using namespace ATOOLS;

namespace {

  //! the provider's (n+np)-leg labels packed: fermion legs as for the
  //! (n+1)-point, photon q in bit nlg+q (set for -1)
  int PackedLabelsN(const std::vector<int> &h, size_t nl, int nlg, size_t np)
  {
    if (h.size() != nl + np) return -1;
    int b(0);
    for (size_t i(0); i < nl; ++i) {
      const int bit(Amplitude::Bit(i));
      if (bit < 0) { if (h[i] != 0) return -1; continue; }
      if (h[i] == 0) return -1;
      if (h[i] < 0) b |= (1 << bit);
    }
    for (size_t q(0); q < np; ++q) {
      if (h[nl + q] == 0) return -1;
      if (h[nl + q] < 0) b |= (1 << (nlg + (int)q));
    }
    return b;
  }

}


bool Acraic_Base::LoopNHelicityToAcraic(const LoopNHelicityFactors &f, int hela,
                                        int helb, std::vector<Complex> &d) const
{
  return LoopNHelicityToAcraic(f, std::vector<int>{hela, helb}, d);
}

bool Acraic_Base::LoopNHelicityToAcraic(const LoopNHelicityFactors &f,
                                        const std::vector<int> &phel,
                                        std::vector<Complex> &d) const
{
  const int nh(Amplitude::NHel()), nlg(Amplitude::s_nlegs);
  const size_t np(phel.size());
  d.assign(nh, Complex(0., 0.));
  if (!f.ok || m_rvmask < 0 || np < 2) return false;
  for (int h : phel) if (h == 0) return false;
  // the (n+1)-point map: its fermion part, and its photon flip on every photon
  const int pf((m_rvmask >> nlg) & 1);
  int mask(m_rvmask & (nh - 1));
  for (size_t q(0); q < np; ++q) mask |= pf << (nlg + (int)q);
  size_t nset(0);
  for (size_t c(0); c < f.hel.size() && c < f.rem.size(); ++c) {
    const int base(PackedLabelsN(f.hel[c], m_flavs.size(), nlg, np));
    if (base < 0) return false;
    const int idx(base ^ mask);
    bool match(true);
    for (size_t q(0); q < np; ++q)
      match = match && ((idx >> (nlg + (int)q)) & 1) == (phel[q] > 0 ? 0 : 1);
    if (!match) continue;
    d[idx & (nh - 1)] = f.rem[c];
    ++nset;
  }
  return nset > 0;
}


double Acraic_Base::LoopNRemainderRho(const std::vector<LoopNHelicityFactors> &pairs,
                                      const std::vector<std::vector<Complex> > &r1)
{
  LoopLadderProbe(pairs, r1);
  if (m_loopbetamax < 2 || pairs.empty()) return 0.;
  const Amplitude &A(m_order == 2 ? m_AmpExpo2 : m_AmpExpo1);
  const Vec4D_Vector &cph(AllPhotonsLab());
  const int nh(Amplitude::NHel());
  auto match = [&](const Vec4D &k) -> long {
    long a(-1); double best(1e30);
    for (size_t c(0); c < cph.size(); ++c) {
      const Vec4D dk(cph[c] - k);
      const double m(Max(Max(dabs(dk[0]), dabs(dk[1])), Max(dabs(dk[2]), dabs(dk[3]))));
      if (m < best) { best = m; a = (long)c; }
    }
    return (a >= 0 && best < 1e-6*Max(k.E(), 1e-30)) ? a : -1;
  };
  double rho(0.);
  for (const LoopNHelicityFactors &f : pairs) {
    if (!f.ok || f.ks.size() != 2) continue;   // pair records only (triples: ladder probe)
    const long a(match(f.k1)), b(match(f.k2));
    if (a < 0 || b < 0 || a == b) { ++m_lpn_nomatch; continue; }
    const PhotonSet S{(int)Min(a, b), (int)Max(a, b)};
    std::map<PhotonSet, Amplitude>::const_iterator iP(m_lpP.find(S)),
      iA(m_lpQa.find(S)), iB(m_lpQb.find(S));
    if (iP == m_lpP.end() || iA == m_lpQa.end() || iB == m_lpQb.end()) continue;
    std::vector<Complex> d;
    if (!LoopNHelicityToAcraic(f, PhotonHelicity(a), PhotonHelicity(b), d)) {
      ++m_lpn_nomap; continue;
    }
    const std::vector<Complex> &ra((size_t)S[0] < r1.size() ? r1[S[0]] : std::vector<Complex>());
    const std::vector<Complex> &rb((size_t)S[1] < r1.size() ? r1[S[1]] : std::vector<Complex>());
    double rpair(0.);
    for (int h = 0; h < nh; ++h) {
      const Complex tD(d[h]*iP->second.m_A[h]);
      const Complex tQ((h < (int)ra.size() ? ra[h] : Complex(0., 0.))*iA->second.m_A[h]
                       + (h < (int)rb.size() ? rb[h] : Complex(0., 0.))*iB->second.m_A[h]);
      const Complex X(tD - tQ);
      rpair += 2.*std::real(conj(A.m_A[h])*X)/4.;
    }
    rho += rpair;
    ++m_lpn_pairs;
  }
  ++m_lpn_ev;
  const double rcr(GetRhoCrude());
  if (rcr > 0.) {
    const double q(rho/rcr);
    m_lpn_sum += q; m_lpn_sum2 += q*q; m_lpn_max = Max(m_lpn_max, std::abs(q));
  }
  return rho;
}


void Acraic_Base::LoopLadderProbe(const std::vector<LoopNHelicityFactors> &pairs,
                                  const std::vector<std::vector<Complex> > &r1)
{
  static const int nev(Settings::GetMainSettings()["ACRAIC"]["LOOP_LADDER_PROBE"].Get<int>());
  static int done(0);
  const size_t n(m_allphotons.size());
  const int nh(Amplitude::NHel());
  if (nev <= 0 || done >= nev || n < 1 || n > 3 || r1.size() < n) return;
  for (size_t j(0); j < n; ++j) if ((int)r1[j].size() != nh) return;
  // delta_T - delta_0 per subset (bit mask of photons), ACRAIC's index: the
  // singles from REAL_VIRTUAL helicity, the larger sets from their records
  std::map<int, std::vector<Complex> > rem;
  for (size_t j(0); j < n; ++j) rem[1 << j] = r1[j];
  if (n >= 2) {
    const Vec4D_Vector &cph(AllPhotonsLab());
    auto same = [](const Vec4D &a, const Vec4D &b) {
      return (a - b).PSpat() + dabs(a[0] - b[0]) < 1e-9*Max(a[0], 1e-30); };
    for (const LoopNHelicityFactors &f : pairs) {
      if (!f.ok || f.ks.size() < 2) continue;
      int T(0);
      std::vector<int> phel;
      bool ok(true);
      for (const Vec4D &k : f.ks) {
        int c(-1);
        for (size_t q(0); q < n && q < cph.size(); ++q) if (same(cph[q], k)) c = (int)q;
        if (c < 0 || (T & (1 << c))) { ok = false; break; }
        T |= 1 << c; phel.push_back(PhotonHelicity(c));
      }
      std::vector<Complex> d;
      if (ok && LoopNHelicityToAcraic(f, phel, d)) rem[T] = d;
    }
    // every subset of two or more photons needs its record
    for (int T(1); T < (1 << n); ++T) if (!rem.count(T)) return;
  }
  // trees on the physical legs, no shifts; soft factors of every leg
  const int all((1 << n) - 1);
  std::map<int, Amplitude> M;
  PropShifts none;
  for (int T(1); T <= all; ++T) {
    PhotonSet S;
    Vec4D_Vector pp(m_pacraic);
    for (size_t j(0); j < n; ++j)
      if (T & (1 << j)) { S.push_back((int)j); pp.push_back(m_allphotons[j]); }
    if (!ComixRealNShifted(pp, S, M[T], none, true)) return;
  }
  std::vector<Complex> s(n);
  for (size_t j(0); j < n; ++j) s[j] = TotalEikonal(m_pacraic, m_allphotons[j], PhotonHelicity(j));
  auto soft = [&](int U) { Complex p(1., 0.);
    for (size_t j(0); j < n; ++j) if (U & (1 << j)) p *= s[j];
    return p; };
  auto bits = [](int T) { int c(0); while (T) { c += T & 1; T >>= 1; } return c; };
  // rho_T, smallest sets first
  std::map<int, std::vector<Complex> > rho;
  for (int size(1); size <= (int)n; ++size)
    for (int T(1); T <= all; ++T) {
      if (bits(T) != size) continue;
      std::vector<Complex> &r(rho[T]);
      r.assign(nh, Complex(0., 0.));
      for (int h(0); h < nh; ++h) r[h] = rem[T][h]*M[T].m_A[h];
      for (int U(1); U < T; ++U) {
        if ((U & T) != U) continue;
        const Complex w(soft(T & ~U));
        for (int h(0); h < nh; ++h) r[h] -= rho[U][h]*w;
      }
    }
  double nM(0.), tgt(0.);
  for (int h(0); h < nh; ++h) {
    nM += std::norm(M[all].m_A[h]);
    tgt += 2.*std::real(conj(M[all].m_A[h])*rem[all][h]*M[all].m_A[h]);
  }
  if (!(nM > 0.)) return;
  ++done;
  std::vector<double> xs;
  const double rs(sqrt(m_s));
  for (size_t j(0); j < n; ++j) xs.push_back(2.*m_allphotons[j][0]/rs);
  std::sort(xs.rbegin(), xs.rend());
  std::ostringstream o;
  o<<std::setprecision(17)<<"@@@ LLCLOS n="<<n<<" x=";
  for (size_t j(0); j < n; ++j) o<<(j ? "," : "")<<xs[j];
  o<<" rem="<<tgt/nM;
  std::vector<Complex> L(nh, Complex(0., 0.));
  for (size_t J(0); J <= n; ++J) {
    if (J > 0)
      for (int T(1); T <= all; ++T) {
        if (bits(T) != (int)J) continue;
        const Complex w(soft(all & ~T));
        for (int h(0); h < nh; ++h) L[h] += rho[T][h]*w;
      }
    double e(0.);
    for (int h(0); h < nh; ++h)
      e += 2.*std::real(conj(M[all].m_A[h])*(L[h] - rem[all][h]*M[all].m_A[h]));
    o<<" d"<<J<<"="<<dabs(e)/nM;
  }
  o<<"\n";
  std::cerr<<o.str();
}


void Acraic_Base::ReportLoopN()
{
  if (m_loopbetamax < 2) return;
  double buf[6] = {(double)m_lpn_ev, (double)m_lpn_pairs, (double)m_lpn_nomatch,
                   (double)m_lpn_nomap, m_lpn_sum, m_lpn_sum2};
  double mx[1] = {m_lpn_max};
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 6, MPI_DOUBLE, MPI_SUM);
    mpi->Allreduce(mx, 1, MPI_DOUBLE, MPI_MAX);
  }
#endif
  if (buf[0] <= 0.) return;
  const double n(buf[0]), mean(buf[4]/n), rms(sqrt(Max(0., buf[5]/n - mean*mean)));
  msg_Out()<<std::setprecision(6)<<"ACRAIC: LOOP_BETA_MAX 2 - "<<(long)buf[0]
           <<" events with photon pairs, "<<(long)buf[1]<<" pair residuals added, "
           <<(long)buf[2]<<" unmatched, "<<(long)buf[3]<<" without a helicity map\n"
           <<"  <rho_loop/rho_crude> per such event: mean "<<mean<<" rms "<<rms
           <<" max |.| "<<mx[0]<<std::endl;
}
