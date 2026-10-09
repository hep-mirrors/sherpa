/*!
  \file Acraic_Beta2.C

  Tree-level residuals beta-hat_m, m >= 2 (ACRAIC: ORDER: 2), from Comix's
  m-photon amplitudes by one subset recursion (ACRAIC: TREE_BETA_MAX).
  Per partition and helicity, for a sorted photon set S,

    beta_S = d_S M_|S|(S) - sum_{T proper subset of S} (prod_{i in S\T} s_i) beta_T,

  beta_{} = m_partborn0, beta_{j} = what the beta_1 pass added, s_i = m_sactu,
  d_S = [all of S on the reducing stage]; the partition sum gets
  (prod_{i not in S} s_i) beta_S. At |S| = 2 this is hep-ph/0006359 eq. (95).
  M_|S| is one Comix amplitude with the other photons as per-stage
  propagator shifts.
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "YFS/NLO/Real_Correction.H"

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"

#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace YFS;
using namespace ATOOLS;

namespace {

  typedef std::vector<int> PhotonSet;

  // photon part of a Comix helicity index (photon q above the fermion legs)
  size_t PhotonHelicityWord(const std::vector<int> &hel, int photonflip)
  {
    size_t word(0);
    for (size_t q(0); q < hel.size(); ++q)
      if (FlipLegs(HelicityBit<int>(hel[q]), photonflip) != 0) word = AddLeg(word, q);
    return word;
  }

  //! Comix's index of fermion configuration ferm with photon word pw
  inline size_t FlatIndexN(size_t ferm, size_t pw, int nlegs)
  { return JoinLegs(ferm, pw, nlegs); }

  //! the number of entries of an amplitude with nlegs fermion and m photon legs
  inline size_t HelicityCount(int nlegs, size_t m)
  { return NumLegSubsets(nlegs + m); }

  void CollectSubsets(const PhotonSet &from, size_t m, size_t start,
                      PhotonSet &current, std::vector<PhotonSet> &out)
  {
    if (current.size() == m) { out.push_back(current); return; }
    for (size_t i(start); i + (m - current.size()) <= from.size(); ++i) {
      current.push_back(from[i]);
      CollectSubsets(from, m, i + 1, current, out);
      current.pop_back();
    }
  }

  //! every m-element subset of the ascending set `from`, in lexicographic order
  std::vector<PhotonSet> SubsetsOfSize(const PhotonSet &from, size_t m)
  {
    std::vector<PhotonSet> out;
    PhotonSet current;
    if (m <= from.size()) CollectSubsets(from, m, 0, current, out);
    return out;
  }

  //! S \ T for ascending S and T
  PhotonSet Without(const PhotonSet &S, const PhotonSet &T)
  {
    PhotonSet d;
    std::set_difference(S.begin(), S.end(), T.begin(), T.end(),
                        std::back_inserter(d));
    return d;
  }

  inline bool Contains(const PhotonSet &S, int i)
  { return std::binary_search(S.begin(), S.end(), i); }

  //! {0, 1, ..., n-1}
  PhotonSet PhotonRange(size_t n)
  {
    PhotonSet a(n);
    for (size_t i(0); i < n; ++i) a[i] = (int)i;
    return a;
  }

}


void Acraic_Base::ResetBeta1Record()
{
  const size_t n(m_allphotons.size());
  if (m_b2M1.size() != n) m_b2M1.assign(n, Amplitude());
  m_b2hasM1.assign(n, 0);
  m_b2on.assign(n, 0);
  m_b2sj.assign(n, Complex(0., 0.));
}


// Why the m-photon real provider cannot serve beta_m (empty if it can).
static std::string TreeProviderRefusal(YFS::Real_Correction *prov,
                                       const Flavour_Vector &flavs, size_t m)
{
  if (prov == NULL || prov->p_proc == NULL)
    return "no "+ToString(m)+"-photon real provider (YFS: NLO_Part must "
           "contain W, and YFS builds the m-photon processes up to "
           "max(YFS: NLO_MAX_PHOTONS, ACRAIC: TREE_BETA_MAX))";
  const Flavour_Vector &pf(prov->p_proc->Flavours());
  const size_t nl(flavs.size());
  if (pf.size() != nl + m) return "the provider's process has the wrong leg count";
  for (size_t i(0); i < nl; ++i)
    if (pf[i] != flavs[i]) return "the provider's process has different legs";
  for (size_t i(nl); i < nl + m; ++i)
    if (pf[i].Kfcode() != kf_photon) return "the provider's extra legs are not photons";
  return "";
}


/*
  Per event: whether beta_m (m >= 2) runs, the eligible photons (not on a
  fixed stage, x = 2E/sqrt(s) > ACRAIC: BETA2_XCUT) and m_treebetacap.
  Refused for Borns with final-state photons (helicity slots unvalidated).
*/
bool Acraic_Base::PrepareBeta2()
{
  const size_t n(m_allphotons.size());
  m_b2elig.assign(n, 0);
  m_treebetacap = 0;
  if (m_treeprobe) m_treeorder.assign(Max(m_treebetamax, 1) + 1, Amplitude());
  if (m_bnterms.size() < (size_t)Max(m_treebetamax, 2) + 1) {
    m_bnterms.resize(Max(m_treebetamax, 2) + 1, 0);
    m_bnamps.resize(Max(m_treebetamax, 2) + 1, 0);
  }
  if (m_order != 2 || n < 2 || m_treebetamax < 2) return false;
  static int state(0);   // 0 unchecked, 1 ok, -1 refused
  static int provcap(2);
  if (state == 0) {
    std::string why(TreeProviderRefusal(RealProvider(2), m_flavs, 2));
    for (size_t i(2); i < m_flavs.size() && why.empty(); ++i)
      if (m_flavs[i].Kfcode() == kf_photon)
        why = "the Born has final-state photons (identical-photon helicity "
              "slots not validated)";
    if (!why.empty()) {
      msg_Error()<<METHOD<<"(): ACRAIC: ORDER 2 requested but beta_2 is OFF for "
                 <<"this run: "<<why<<". The ACRAIC column is O(alpha^1)."<<std::endl;
      state = -1;
    }
    else {
      state = 1;
      msg_Info()<<"ACRAIC: beta_2 from "<<RealProvider(2)->p_proc->Name()<<std::endl;
      for (int m(3); m <= m_treebetamax; ++m) {
        const std::string whym(TreeProviderRefusal(RealProvider(m), m_flavs, m));
        if (!whym.empty()) {
          msg_Error()<<METHOD<<"(): ACRAIC: TREE_BETA_MAX "<<m_treebetamax
                     <<" requested, but "<<whym<<"; the tree residuals stop at "
                     <<"beta_"<<(m - 1)<<"."<<std::endl;
          break;
        }
        provcap = m;
        msg_Info()<<"ACRAIC: beta_"<<m<<" from "<<RealProvider(m)->p_proc->Name()
                  <<std::endl;
      }
    }
  }
  if (state < 0) return false;
  static const double xcut(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                           ["BETA2_XCUT"].Get<double>());
  const double rs(m_s > 0. ? sqrt(m_s) : 0.);
  if (!(rs > 0.)) return false;
  size_t ne(0);
  for (size_t i(0); i < n; ++i) {
    // a collapsed (soft) photon has no beta_2; without FSR every photon is
    // "fixed" only because there is one stage to put it on (PartitionStart)
    // - nu nu, ISR only - and before 2026-10-08 that silently left those
    // processes with no beta_2 at ORDER 2
    if (HasFSR() && m_fixedstage.size() == n && m_fixedstage[i]) continue;
    if (2.*m_allphotons[i][0]/rs > xcut) { m_b2elig[i] = 1; ++ne; }
  }
  if (ne < 2) return false;
  m_treebetacap = Min(Min(provcap, m_treebetamax), (int)ne);
  return true;
}


bool Acraic_Base::Beta2PairEligible(int j, int l) const
{
  return j >= 0 && l >= 0 && j != l && j < (int)m_b2elig.size()
    && l < (int)m_b2elig.size() && m_b2elig[j] && m_b2elig[l];
}


Acraic_Base::PhotonSet Acraic_Base::TreeBetaPhotons() const
{
  PhotonSet e;
  for (size_t i(0); i < m_b2elig.size(); ++i)
    if (m_b2elig[i]) e.push_back((int)i);
  return e;
}


// Comix's m-photon amplitude in the partition-sum alignment; hel[q] is the
// helicity of photon leg q.
bool Acraic_Base::ComixRealNHel(const Vec4D_Vector &pp, const std::vector<int> &hel,
                              Amplitude &M, const PropShifts &shifts,
                              bool noshift)
{
  const size_t m(hel.size());
  YFS::Real_Correction *prov(RealProvider(m));
  if (prov == NULL) return false;
  const int nlg(Amplitude::s_nlegs);
  if (pp.size() != m_flavs.size() + m) return false;
  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (noshift ? prov->ComixAmplitudes(pp, -1.)
             : prov->ComixAmplitudesShifts(pp, shifts));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if (sa.size() != HelicityCount(nlg, m)) return false;
  const int pf(ComixPhotonFlip(nlg));
  const size_t pw(PhotonHelicityWord(hel, pf));
  const Complex cpl(ComixPhotonCoupling(m), 0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f)
    M.m_A[f] = m_cxbalign.m_A[f] * cpl
      * sa[FlatIndexN((size_t)FlipLegs(f, ComixLegFlips(nlg)), pw, nlg)];
  return true;
}

bool Acraic_Base::ComixRealNShifted(const Vec4D_Vector &pp, const PhotonSet &S,
                                  Amplitude &M, const PropShifts &shifts,
                                  bool noshift)
{
  std::vector<int> hel;
  for (size_t q(0); q < S.size(); ++q) hel.push_back(m_PhoHel[S[q]]);
  return ComixRealNHel(pp, hel, M, shifts, noshift);
}

bool Acraic_Base::ComixReal2Shifted(const Vec4D_Vector &pp, int hj, int hl,
                                  Amplitude &M2, const PropShifts &shifts,
                                  bool noshift)
{
  return ComixRealNHel(pp, std::vector<int>{hj, hl}, M2, shifts, noshift);
}


Vec4D Acraic_Base::StagePhotonSumExcept(int stage, const PhotonSet &S) const
{
  Vec4D K;
  if (m_stage.size() != m_allphotons.size()) return K;
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if (!Contains(S, (int)i) && m_stage[i] == stage) K += m_allphotons[i];
  return K;
}


Acraic_Base::PropShifts Acraic_Base::SubsetStageShifts(const PhotonSet &S) const
{
  PropShifts sh;
  if (m_stage.size() != m_allphotons.size()) return sh;
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    const size_t mask(StageShiftMask((int)g));   // as StageShifts
    if (mask == 0) continue;
    sh.push_back(std::make_pair(mask, StagePhotonSumExcept((int)g, S)));
  }
  AddResonanceEntries(sh);   // ACRAIC: RESONANCE_EXPONENTIATION (else none)
  return sh;
}

Acraic_Base::PropShifts Acraic_Base::StageShifts2(int j, int l) const
{
  return SubsetStageShifts(PhotonSet{Min(j, l), Max(j, l)});
}


bool Acraic_Base::SubsetPartitionLegs(const PhotonSet &S, Vec4D_Vector &pb) const
{
  const int n((int)m_allphotons.size());
  if (m_stage.size() != m_allphotons.size()) return false;
  for (size_t q(0); q < S.size(); ++q)
    if (S[q] < 0 || S[q] >= n) return false;
  Vec4D X(m_pacraic[0] + m_pacraic[1]);
  for (int i(0); i < n; ++i)
    if (!Contains(S, i) && m_stagereduces[m_stage[i]]) X -= m_allphotons[i];
  Vec4D Y(X);
  for (size_t q(0); q < S.size(); ++q) Y -= m_allphotons[S[q]];
  for (size_t i(2); i < m_pacraic.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pacraic[i];
  return LegsAt(X, Y, pb);
}


void Acraic_Base::AddSubsetExchangeLineShifts(const PhotonSet &S, PropShifts &sh) const
{
  if (!ExchangeLineShiftsOn() || m_pacraic.size() < 4
      || m_flavs.size() != m_pacraic.size()) return;
  Vec4D_Vector pb;
  if (!SubsetPartitionLegs(S, pb) || pb.size() != m_pacraic.size()) return;
  for (size_t i(0); i < m_pacraic.size(); ++i) {
    const Vec4D d((i < 2 ? -1. : 1.)*(pb[i] - m_pacraic[i]));
    sh.push_back(std::make_pair(LegUnion(LegBit(i),
                                         PHASIC::Process_Base::s_propshiftleg), d));
  }
}


Complex Acraic_Base::RecordedBeta1(int j, int f) const
{
  const Complex B(m_partborn0.m_A[f]);
  return m_b2on[j] ? (m_b2hasM1[j] ? m_b2M1[j].m_A[f] : Complex(0., 0.))
                     - m_b2sj[j]*B : Complex(0., 0.);
}


// as a product, never sProd/prod, which has poles at the zeros of the s_i
Complex Acraic_Base::SoftProduct(const PhotonSet &S) const
{
  Complex w(1., 0.);
  for (size_t q(0); q < S.size(); ++q) w *= m_sactu[S[q]];
  return w;
}


/*
  One set S in the current partition: W_S [ d_S M_|S| - sum_T (prod_{S\T} s) beta_T ],
  subsets largest first, which at |S| = 2 is eq. (95) in the order
  M_2 - s_l beta_1(k_j) - s_j beta_1(k_l) - s_j s_l B.
*/
bool Acraic_Base::ComixInfraredSubtracted_N_0(const PhotonSet &S)
{
  const size_t n(m_allphotons.size()), m(S.size());
  if (!m_cxbalignok || m < 2 || m_stage.size() != n
      || m_sactu.size() != n || m_b2on.size() != n) return false;
  std::vector<const Amplitude*> lower;
  std::vector<Complex> coef;
  for (size_t t(m); t-- > 0; ) {
    const std::vector<PhotonSet> sub(SubsetsOfSize(S, t));
    for (size_t a(0); a < sub.size(); ++a) {
      std::map<PhotonSet, Amplitude>::const_iterator it(m_treebeta.find(sub[a]));
      if (it == m_treebeta.end()) { ++m_b2fail; return false; }
      lower.push_back(&it->second);
      coef.push_back(SoftProduct(Without(S, sub[a])));
    }
  }
  bool allreducing(true);
  for (size_t q(0); q < m; ++q)
    if (!m_stagereduces[m_stage[S[q]]]) allreducing = false;
  // YFS: IFI false: every partition carries its placement's share of M_|S|
  const bool addm(!m_ifi || allreducing);
  Amplitude M;
  if (addm) {
    static const beta1legs::code legsmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                              ["BETA1_LEGS"].Get<beta1legs::code>());
    const bool rebuild(legsmode == beta1legs::balanced || m_redborn);
    Vec4D_Vector pb;
    PropShifts shifts;
    if (rebuild) { if (!SubsetPartitionLegs(S, pb)) { ++m_b2fail; return false; } }
    else {
      pb = m_pacraic;
      shifts = SubsetStageShifts(S);
      AddSubsetExchangeLineShifts(S, shifts);
      if (!m_ifi || m_isronly) AddPlacementVetoes(S, shifts);
    }
    Vec4D_Vector pp(pb);
    for (size_t q(0); q < m; ++q) pp.push_back(m_allphotons[S[q]]);
    if (!ComixRealNShifted(pp, S, M, shifts, rebuild)) {
      ++m_b2fail;
      static bool warned(false);
      if (!warned) {
        warned = true;
        msg_Error()<<METHOD<<"(): Comix "<<m<<"-photon amplitude unavailable; "
                   <<"this set's beta_"<<m<<" is dropped. Reported once."<<std::endl;
      }
      return false;
    }
    ++m_bnamps[m];
  }
  ++m_bnterms[m];
  const Complex w(SoftProduct(Without(PhotonRange(n), S)));
  Amplitude &beta(m_treebeta[S]);
  const int nh(Amplitude::NHel());
  double nv(0.), nm(0.), nss(0.);
  for (int f = 0; f < nh; ++f) {
    Complex v(addm ? M.m_A[f] : Complex(0., 0.));
    for (size_t a(0); a < lower.size(); ++a) v -= coef[a]*lower[a]->m_A[f];
    beta.m_A[f] = v;
    const Complex add(w*v);
    m_AmpBeta2.m_A[f] += add;
    if (m_treeprobe && m < m_treeorder.size()) m_treeorder[m].m_A[f] += add;
    nv += std::norm(add);
    if (addm) nm += std::norm(w*M.m_A[f]);
    nss += std::norm(w*coef.back()*lower.back()->m_A[f]);   // W prod_S s B
  }
  static const int tr(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                      ["BETA2_TRACE"].Get<int>());
  static int ntr(0);
  if (tr && ntr < tr*64) {
    ++ntr;
    std::string st, set;
    for (size_t i(0); i < m_stage.size(); ++i) st += ToString(m_stage[i]);
    for (size_t q(0); q < m; ++q) set += (q ? "," : "") + ToString(S[q]);
    std::ostringstream o;
    o<<std::setprecision(6)<<(m == 2 ? "@@@ B2TRACE" : "@@@ BNTRACE")
     <<" wp="<<st<<" S="<<set<<" m="<<m<<" all_reducing="<<allreducing
     <<" |W prod_S s B|="<<sqrt(nss)<<" |W M|="<<sqrt(nm)
     <<" |W beta|="<<sqrt(nv)<<"\n";
    std::cerr<<o.str();
  }
  return true;
}


// Probe ladder rungs m = 0, 1: W_{} beta_{} and W_{j} beta_{j}.
void Acraic_Base::AddProbeLowOrders()
{
  if (m_treeorder.size() < 2) return;
  const size_t n(m_allphotons.size());
  const PhotonSet all(PhotonRange(n));
  const int nh(Amplitude::NHel());
  const Complex w0(SoftProduct(all));
  for (int f = 0; f < nh; ++f) m_treeorder[0].m_A[f] += w0*m_partborn0.m_A[f];
  for (size_t j(0); j < n; ++j) {
    const Amplitude &b(m_treebeta[PhotonSet{(int)j}]);
    const Complex w(SoftProduct(Without(all, PhotonSet{(int)j})));
    for (int f = 0; f < nh; ++f) m_treeorder[1].m_A[f] += w*b.m_A[f];
  }
}


// Every beta_S of the current partition, smallest sets first.
void Acraic_Base::AddTreeBetas(bool dobetan)
{
  const size_t n(m_allphotons.size());
  if (m_b2on.size() != n) return;
  m_treebeta.clear();
  m_treebeta[PhotonSet()] = m_partborn0;
  const int nh(Amplitude::NHel());
  for (size_t j(0); j < n; ++j) {
    Amplitude &b(m_treebeta[PhotonSet{(int)j}]);
    for (int f = 0; f < nh; ++f) b.m_A[f] = RecordedBeta1((int)j, f);
  }
  if (m_treeprobe) AddProbeLowOrders();
  if (!dobetan) return;
  const PhotonSet elig(TreeBetaPhotons());
  for (int m(2); m <= m_treebetacap; ++m) {
    const std::vector<PhotonSet> sets(SubsetsOfSize(elig, (size_t)m));
    for (size_t a(0); a < sets.size(); ++a) ComixInfraredSubtracted_N_0(sets[a]);
  }
}


// ACRAIC: TREE_BETA_PROBE: A_N = sum_{m <= N} m_treeorder[m] against Comix's
// unshifted n-photon tree; exact when every photon is eligible and cap >= n.
void Acraic_Base::TreeBetaClosure()
{
  static const int nev(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["TREE_BETA_PROBE"].Get<int>());
  static int done(0);
  const size_t n(m_allphotons.size());
  if (!m_treeprobe || done >= nev || n < 1 || !m_cxbalignok
      || m_treeorder.size() < 2) return;
  ++done;
  Vec4D_Vector pp(m_pacraic);
  for (size_t j(0); j < n; ++j) pp.push_back(m_allphotons[j]);
  Amplitude Mex;
  PropShifts none;
  const bool haveM(RealProvider(n) != NULL
                   && ComixRealNShifted(pp, PhotonRange(n), Mex, none, true));
  const int nh(Amplitude::NHel());
  double nM(0.), maxM(0.);
  for (int f = 0; f < nh; ++f) {
    nM += std::norm(Mex.m_A[f]);
    maxM = Max(maxM, std::abs(Mex.m_A[f]));
  }
  size_t ne(0), nfixed(0);
  for (size_t j(0); j < n; ++j) {
    if (j < m_b2elig.size() && m_b2elig[j]) ++ne;
    if (m_fixedstage.size() == n && m_fixedstage[j]) ++nfixed;
  }
  // per photon x and rounding floor eps E_l^2/(p_l.k), max over charged legs
  std::vector<std::pair<double, double> > xf;
  const double rs(sqrt(m_s));
  for (size_t j(0); j < n; ++j) {
    const Vec4D &k(m_allphotons[j]);
    double floor(0.);
    for (size_t i(0); i < m_pacraic.size(); ++i) {
      if (m_flavs[i].Charge() == 0.) continue;
      const double pk(m_pacraic[i]*k);
      if (pk > 0.) floor = Max(floor, 1.1e-16*sqr(m_pacraic[i][0])/pk);
    }
    xf.push_back(std::make_pair(2.*k[0]/rs, floor));
  }
  std::sort(xf.rbegin(), xf.rend());
  std::vector<double> xs;
  for (size_t j(0); j < n; ++j) xs.push_back(xf[j].first);
  std::ostringstream o;
  o<<std::setprecision(17)<<"@@@ BNCLOS n="<<n<<" nI="<<m_isrphotons.size()
   <<" nF="<<m_fsrphotons.size()<<" elig="<<(n == 1 ? 1 : ne)
   <<" fixed="<<nfixed<<" cap="<<(n == 1 ? 1 : m_treebetacap)<<" x=";
  for (size_t j(0); j < n; ++j) o<<(j ? "," : "")<<xs[j];
  o<<" floor=";
  for (size_t j(0); j < n; ++j) o<<(j ? "," : "")<<std::setprecision(3)<<xf[j].second;
  o<<std::setprecision(17)<<" M2="<<(haveM ? nM : -1.);
  Amplitude A;
  const size_t top(Min(m_treeorder.size() - 1, n));
  double nA1(0.), nd1(0.);
  for (size_t N(0); N <= top; ++N) {
    double nA(0.), worst(0.);
    for (int f = 0; f < nh; ++f) {
      A.m_A[f] += m_treeorder[N].m_A[f];
      nA += std::norm(A.m_A[f]);
      // live helicities only: the mass-suppressed ones have |M| ~ 0
      if (haveM && std::abs(Mex.m_A[f]) > 1e-6*maxM)
        worst = Max(worst, std::abs(A.m_A[f] - Mex.m_A[f])/std::abs(Mex.m_A[f]));
      if (N == 1) { nA1 += std::norm(A.m_A[f]);
        nd1 += std::norm(A.m_A[f] - m_AmpExpo1.m_A[f]); }
    }
    o<<" A"<<N<<"="<<nA<<":"<<(haveM ? worst : -1.);
  }
  o<<" A1chk="<<(nA1 > 0. ? sqrt(nd1/nA1) : -1.)<<"\n";
  std::cerr<<o.str();
}


// n = 2 closure: sum_wp [s s B + s beta_1 + s beta_1 + beta_2] = M_2.
void Acraic_Base::Beta2Closure()
{
  if (m_order != 2 || m_allphotons.size() != 2 || !Beta2PairEligible(0, 1)
      || !m_cxbalignok) return;
  Vec4D_Vector pp(m_pacraic);
  pp.push_back(m_allphotons[0]);
  pp.push_back(m_allphotons[1]);
  Amplitude Mex;
  PropShifts none;
  if (!ComixReal2Shifted(pp, m_PhoHel[0], m_PhoHel[1], Mex, none, true)) return;
  double nd(0.), ne(0.), n2(0.), n1(0.), nb(0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) {
    nd += std::norm(m_AmpExpo2.m_A[f] - Mex.m_A[f]);
    ne += std::norm(Mex.m_A[f]);
    n2 += std::norm(m_AmpExpo2.m_A[f]);
    n1 += std::norm(m_AmpExpo1.m_A[f]);
    nb += std::norm(m_AmpBeta2.m_A[f]);
  }
  const double rs(sqrt(m_s));
  const size_t ni(m_isrphotons.size());
  std::ostringstream o;
  o<<std::setprecision(10)<<"@@@ B2CLOS cls="<<(ni == 2 ? "II" : ni == 1 ? "IF" : "FF")
   <<" x0="<<2.*m_allphotons[0][0]/rs<<" x1="<<2.*m_allphotons[1][0]/rs
   <<" h="<<m_PhoHel[0]<<","<<m_PhoHel[1]
   <<" rel="<<(ne > 0. ? sqrt(nd/ne) : -1.)
   <<" rho2="<<n2/4.<<" |M2|^2/4="<<ne/4.<<" rho1="<<n1/4.
   <<" |B2|^2/4="<<nb/4.<<" rho0="<<m_result0<<" nparts="<<m_nparts<<"\n";
  std::cerr<<o.str();
}


// ACRAIC: BETA2_SOFT_TEST: scale the photons by lambda on balanced legs;
// |beta_2|/|M_2| must fall linearly. s-channel 2 -> 2 only.
void Acraic_Base::Beta2SoftLimit()
{
  static const int nev(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BETA2_SOFT_TEST"].Get<int>());
  static int done(0);
  if (m_order != 2 || done >= nev || m_allphotons.size() != 2
      || !Beta2PairEligible(0, 1) || !m_cxbalignok) return;
  ++done;
  const std::vector<int> stsave(m_stage);
  const Vec4D_Vector phsave(m_allphotons);
  const Vec4D P(m_pacraic[0] + m_pacraic[1]);
  const int nh(Amplitude::NHel());
  const double rn(RealNorm());
  const double rs(sqrt(m_s));
  static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
  const size_t ni(m_isrphotons.size());
  const char *cls(ni == 2 ? "II" : ni == 1 ? "IF" : "FF");
  for (int mode(0); mode < 3; ++mode) {
    std::ostringstream o;
    o<<std::setprecision(6)<<"@@@ B2SOFT cls="<<cls<<" mode="
     <<(mode == 0 ? "k1soft" : mode == 1 ? "k0soft" : "both")
     <<" x0="<<2.*phsave[0][0]/rs<<" x1="<<2.*phsave[1][0]/rs
     <<" pfmode="<<pfmode;
    // each photon's angle to its nearest charged leg, also in units of m/E
    for (int j(0); j < 2; ++j) {
      double best(1e30), bu(1e30);
      for (size_t i(0); i < m_pacraic.size(); ++i) {
        if (m_flavs[i].Charge() == 0.) continue;
        const double c(Vec3D(phsave[j])*Vec3D(m_pacraic[i])
                       /(Vec3D(phsave[j]).Abs()*Vec3D(m_pacraic[i]).Abs()));
        const double th(acos(Max(-1., Min(1., c))));
        const double u(m_flavs[i].Mass() > 0. ? th/(m_flavs[i].Mass()/m_pacraic[i][0]) : th);
        if (th < best) { best = th; bu = u; }
      }
      o<<" th"<<j<<"="<<best<<"("<<bu<<"m/E)";
    }
    for (int e(0); e <= 5; ++e) {
      const double lam(pow(10., -(double)e));
      const double l0(mode == 0 ? 1. : lam), l1(mode == 1 ? 1. : lam);
      const Vec4D k0(l0*phsave[0]), k1(l1*phsave[1]);
      Vec4D Y(P - k0 - k1);
      for (size_t i(2); i < m_pacraic.size(); ++i)
        if (i != m_if1 && i != m_if2) Y -= m_pacraic[i];
      Vec4D_Vector pl;
      if (!LegsAt(P, Y, pl)) continue;
      m_allphotons[0] = k0; m_allphotons[1] = k1;
      Amplitude B2, M2s, SSB, B1s, SB;
      bool ok(true);
      for (int w0(0); w0 < m_nstages && ok; ++w0)
        for (int w1(0); w1 < m_nstages && ok; ++w1) {
          m_stage[0] = w0; m_stage[1] = w1;
          const bool d0(m_stagereduces[w0] != 0), d1(m_stagereduces[w1] != 0);
          Amplitude C, M10, M11, M2;
          ok = ComixBornShifted(pl, C, StageShifts(-1));
          const Complex s0(StageEikonal(w0, pl, k0, m_PhoHel[0], 0));
          const Complex s1(StageEikonal(w1, pl, k1, m_PhoHel[1], 1));
          if (ok && d0) { Vec4D_Vector pp(pl); pp.push_back(k0);
            ok = ComixRealShifted(pp, m_PhoHel[0], M10, StageShifts(0)); }
          if (ok && d1) { Vec4D_Vector pp(pl); pp.push_back(k1);
            ok = ComixRealShifted(pp, m_PhoHel[1], M11, StageShifts(1)); }
          if (ok && d0 && d1) { Vec4D_Vector pp(pl); pp.push_back(k0); pp.push_back(k1);
            ok = ComixReal2Shifted(pp, m_PhoHel[0], m_PhoHel[1], M2, StageShifts2(0, 1)); }
          if (!ok) break;
          for (int f = 0; f < nh; ++f) {
            const Complex B(m_cxbalign.m_A[f]*C.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
            const Complex b10((d0 ? m_cxbalign.m_A[f]*M10.m_A[f]/rn : Complex(0., 0.)) - s0*B);
            const Complex b11((d1 ? m_cxbalign.m_A[f]*M11.m_A[f]/rn : Complex(0., 0.)) - s1*B);
            const Complex m2(d0 && d1 ? M2.m_A[f] : Complex(0., 0.));
            B2.m_A[f]  += m2 - s1*b10 - s0*b11 - s0*s1*B;
            M2s.m_A[f] += m2;
            SSB.m_A[f] += s0*s1*B;
            // beta_1/(s_f B) must level off at the LBK constant
            const Complex bs(mode == 1 ? b10 : b11), sf(mode == 1 ? s1 : s0);
            B1s.m_A[f] += sf*bs;
            SB.m_A[f]  += sf*B;
          }
        }
      if (!ok) { o<<" lam="<<lam<<":fail"; continue; }
      double nb(0.), nm(0.), nss(0.), nb1(0.), nsb(0.);
      for (int f = 0; f < nh; ++f) { nb += std::norm(B2.m_A[f]); nm += std::norm(M2s.m_A[f]);
        nss += std::norm(SSB.m_A[f]); nb1 += std::norm(B1s.m_A[f]); nsb += std::norm(SB.m_A[f]); }
      // lam:|beta_2|/|M_2|:|beta_2|/|s s B|:|beta_1|/|s B| (last single-soft only)
      o<<" lam="<<lam<<":"<<(nm > 0. ? sqrt(nb/nm) : -1.)
       <<":"<<(nss > 0. ? sqrt(nb/nss) : -1.);
      if (mode < 2) o<<":"<<(nsb > 0. ? sqrt(nb1/nsb) : -1.);
    }
    o<<"\n";
    std::cerr<<o.str();
  }
  m_stage = stsave;
  m_allphotons = phsave;
}


void Acraic_Base::ReportBeta2() const
{
  if (m_order != 2) return;
  msg_Info()<<"ACRAIC ORDER 2: tree residuals on "<<m_b2events<<" events, "
            <<m_b2fail<<" failures."<<std::endl;
  for (size_t m(2); m < m_bnterms.size(); ++m)
    if (m_bnterms[m] > 0)
      msg_Info()<<"  beta_"<<m<<": "<<m_bnterms[m]<<" set-partition terms, "
                <<m_bnamps[m]<<" "<<m<<"-photon amplitudes."<<std::endl;
}
