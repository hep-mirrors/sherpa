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
#include <set>
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
  m_b2onept.assign(n, 0);
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
  const size_t nl(flavs.size()), nin(prov->p_proc->NIn());
  if (pf.size() != nl + m) return "the provider's process has the wrong leg count";
  for (size_t i(0); i < nin; ++i)
    if (pf[i] != flavs[i]) return "the provider's process has different legs";
  /*
    The final legs in any order: Sherpa sorts them (a photon goes before
    massless quarks, e- e+ P P d db), and Real_Correction maps the Born legs
    plus photons onto them (m_perm).
  */
  std::multiset<long> want, have;
  for (size_t i(nin); i < nl; ++i) want.insert((long)flavs[i]);
  for (size_t j(0); j < m; ++j) want.insert((long)Flavour(kf_photon));
  for (size_t i(nin); i < pf.size(); ++i) have.insert((long)pf[i]);
  if (want != have) return "the provider's process has different legs";
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
  m_lpP.clear(); m_lpQa.clear(); m_lpQb.clear();
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
    /*
      ACRAIC: BETA2_BORN_PHOTONS (default false until validated): beta_m
      for Borns with final-state photons. The m-photon amplitudes are the
      full Comix ones (symmetric in all photons, no symmetry factor), and
      the Born-pair multichannel weight (YFS: REAL_BORN_PHOTON_MULTICHANNEL,
      rho_crude x G over all Born-pair choices) covers any number of ISR
      photons, so the open question is only the helicity-slot alignment of
      the extra photons next to the Born ones.
    */
    /*
      With TCHANNEL_MULTIPHOTON factorised the pairs are not this recursion
      but AddFactorisedBeta2, at each pair's own scaled point, whose total at
      two photons is M_2 by construction: allowed without the switch.
    */
    static const bool bornphotons(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                  ["BETA2_BORN_PHOTONS"].SetDefault(false).Get<bool>()
                                  || TchannelMultiphotonMode() == tchmultiphoton::factorised);
    for (size_t i(2); i < m_flavs.size() && why.empty() && !bornphotons; ++i)
      if (m_flavs[i].Kfcode() == kf_photon)
        why = "the Born has final-state photons (identical-photon helicity "
              "slots not validated; ACRAIC: BETA2_BORN_PHOTONS true to try)";
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
  /*
    A collapsed (soft) photon has no beta_2. With a single radiating stage
    every photon is "fixed" only because there is one stage to put it on
    (PartitionStart), and its total eikonal is that stage's: nu nu, MODE ISR,
    and MODE FSR with one final stage. Those fixed photons stay eligible.
  */
  int nradiating(0);
  for (int g(0); g < m_nstages && g < (int)m_stagelegs.size(); ++g)
    for (size_t l(0); l < m_stagelegs[g].size(); ++l)
      if (m_stagelegs[g][l].w != 0.) { ++nradiating; break; }
  size_t ne(0);
  for (size_t i(0); i < n; ++i) {
    if (nradiating > 1 && m_fixedstage.size() == n && m_fixedstage[i]) continue;
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
  /*
    M_|S| goes on the one partition per assignment of the others with every
    enumerated photon of S on a reducing stage. A fixed photon is never
    enumerated, so it counts as placed (as for M_1 in
    ComixInfraredSubtracted_1_0): MODE FSR fixes every photon on a final stage.
  */
  bool placed(true);
  for (size_t q(0); q < m; ++q) {
    const bool fixed(m_fixedstage.size() == m_stage.size() && m_fixedstage[S[q]]);
    if (!fixed && !m_stagereduces[m_stage[S[q]]]) placed = false;
  }
  // YFS: IFI false: every partition carries its placement's share of M_|S|
  const bool addm(!m_ifi || placed);
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
      if (!m_ifi || SingleSided()) AddPlacementVetoes(S, shifts);
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
  const int tnpclass(m_tnp ? TnpSetClass(S) : tnp_b2m);
  // ACRAIC: LOOP_BETA_MAX 2: the pair's partition sums for its one-loop
  // residual (Acraic_LoopN.C)
  const bool loopacc(m_loopbetamax >= 2 && m == 2);
  Amplitude *lpP(loopacc ? &m_lpP[S] : NULL), *lpQa(loopacc ? &m_lpQa[S] : NULL),
            *lpQb(loopacc ? &m_lpQb[S] : NULL);
  for (int f = 0; f < nh; ++f) {
    if (loopacc) {
      const Complex Mj(m_b2hasM1[S[0]] ? m_b2M1[S[0]].m_A[f] : Complex(0., 0.));
      const Complex Ml(m_b2hasM1[S[1]] ? m_b2M1[S[1]].m_A[f] : Complex(0., 0.));
      if (addm) lpP->m_A[f] += w*M.m_A[f];
      lpQa->m_A[f] += w*m_sactu[S[1]]*Mj;
      lpQb->m_A[f] += w*m_sactu[S[0]]*Ml;
    }
    Complex v(addm ? M.m_A[f] : Complex(0., 0.));
    for (size_t a(0); a < lower.size(); ++a) v -= coef[a]*lower[a]->m_A[f];
    beta.m_A[f] = v;
    const Complex add(w*v);
    m_AmpBeta2.m_A[f] += add;
    TnpAdd(tnpclass, f, add);
    if (m_treeprobe && m < m_treeorder.size()) m_treeorder[m].m_A[f] += add;
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
  if (TchannelMultiphotonActive()
      && TchannelMultiphotonMode() == tchmultiphoton::factorised) {
    AddFactorisedBeta2();
    return;
  }
  const PhotonSet elig(TreeBetaPhotons());
  for (int m(2); m <= m_treebetacap; ++m) {
    const std::vector<PhotonSet> sets(SubsetsOfSize(elig, (size_t)m));
    for (size_t a(0); a < sets.size(); ++a) ComixInfraredSubtracted_N_0(sets[a]);
  }
}


/*
  OnePhotonScaledLegs for a set of photons: the event with only the photons
  of S, all scaled by x = sqrt(Q^2/(P - K_S)^2) in the P rest frame so the
  final state keeps Q^2, and carried back by the generator's boost Poincare(Q).
  At |S| = n it is the event itself; as a photon of S goes soft it becomes the
  point of S without it.
*/
bool Acraic_Base::SubsetScaledLegs(const PhotonSet &S, Vec4D_Vector &legs,
                                   Vec4D_Vector &kpt) const
{
  if (S.empty() || m_pacraic.size() != m_flavs.size() || m_pacraic.size() < 4)
    return false;
  for (size_t q(0); q < S.size(); ++q)
    if (S[q] < 0 || S[q] >= (int)m_allphotons.size()) return false;
  const Vec4D P(m_pacraic[0] + m_pacraic[1]);
  Vec4D Q, K;
  for (size_t i(2); i < m_pacraic.size(); ++i) Q += m_pacraic[i];
  for (size_t q(0); q < S.size(); ++q) K += m_allphotons[S[q]];
  const Vec4D R(P - K);
  const double P2(P.Abs2()), Q2(Q.Abs2()), R2(R.Abs2());
  if (!(P2 > 0.) || !(Q2 > 0.) || !(R2 > 0.) || !(Q[0] > 0.)) return false;
  const double x(Min(1., sqrt(Q2/R2)));
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  const double sx(x*x*P2);
  const double kal(sqr(sx - ma*ma - mb*mb) - 4.*ma*ma*mb*mb);
  if (!(sx > sqr(ma + mb)) || !(kal > 0.)) return false;
  Poincare toP(P), toQ(Q);
  kpt.assign(S.size(), Vec4D());
  Vec4D Kp;
  for (size_t q(0); q < S.size(); ++q) {
    Vec4D k(m_allphotons[S[q]]);
    toP.Boost(k);
    k *= x;
    kpt[q] = k;
    Kp += k;
  }
  const double pz(sqrt(kal)/(2.*sqrt(sx)));
  const Vec4D pa(sqrt(pz*pz + ma*ma), 0., 0.,  pz);
  const Vec4D pb(sqrt(pz*pz + mb*mb), 0., 0., -pz);
  const Vec4D Qp(pa + pb - Kp);
  if (!(Qp[0] > 0.) || !(Qp.Abs2() > 0.)) return false;
  Poincare toQp(Qp);
  legs = m_pacraic;
  legs[0] = pa; legs[1] = pb;
  for (size_t i(2); i < m_pacraic.size(); ++i) {
    Vec4D p(m_pacraic[i]);
    toQ.Boost(p);
    toQp.BoostBack(p);
    legs[i] = p;
  }
  for (Vec4D &l : legs) toQ.BoostBack(l);
  for (Vec4D &k : kpt) toQ.BoostBack(k);
  Vec4D bal(legs[0] + legs[1]);
  for (size_t i(2); i < legs.size(); ++i) bal -= legs[i];
  for (const Vec4D &k : kpt) bal -= k;
  const double scale(Max(legs[0][0] + legs[1][0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-10*scale) return false;
  return true;
}


/*
  ACRAIC: TCHANNEL_MULTIPHOTON factorised at ORDER 2. In this mode photon j's
  beta_1 is s_j (u_j - B_0), u_j = M_1/s_j at its one-photon point, and
  AddFactorisedRemainder makes the beta_1 sum the factorised
      T = W [ B_0 prod_j c_j + sum_j (u_j - c_j B_0) prod_{l != j} c_l ],
  c_j = <B_0|u_j>/<B_0|B_0>, W = prod_i s_i (m_Sprod). T already contains a
  second-order term, so the recursion's beta_2 = M_2 - ... would count it
  twice. The pair term is instead taken against T, at the pair's own point:
      Delta_jl = W [ u_jl - u_j c_l - u_l c_j + c_j c_l B_0 ] prod_{m != j,l} c_m,
  u_jl = M_2/(s_j s_l) at SubsetScaledLegs({j, l}). With two photons the
  total is W u_jl = M_2 (the point is the event); as either photon goes soft,
  u_jl -> u_j or u_l and Delta_jl -> 0. Pairs only (beta_3 and higher are not
  built in this mode).
*/
void Acraic_Base::AddFactorisedBeta2()
{
  const size_t n(m_allphotons.size());
  if (n < 2 || m_b2on.size() != n || m_b2onept.size() != n
      || m_sactu.size() != n || !m_cxbalignok) return;
  const int nh(Amplitude::NHel());
  const Amplitude &B0(m_partborn0);
  double nb(0.);
  for (int f(0); f < nh; ++f) nb += std::norm(B0.m_A[f]);
  if (!(nb > 0.) || IsBad(nb)) return;
  // u_j and c_j; a photon without a one-photon-point beta_1 has u_j = B_0, c_j = 1
  std::vector<Amplitude> u(n, B0);
  std::vector<Complex> c(n, Complex(1., 0.));
  std::vector<char> ok(n, 0);
  for (size_t j(0); j < n; ++j) {
    if (!m_b2on[j] || !m_b2hasM1[j] || !m_b2onept[j] || !(std::abs(m_b2sj[j]) > 0.))
      continue;
    Complex proj(0., 0.);
    for (int f(0); f < nh; ++f) {
      u[j].m_A[f] = m_b2M1[j].m_A[f]/m_b2sj[j];
      proj += std::conj(B0.m_A[f])*u[j].m_A[f];
    }
    c[j] = proj/nb;
    ok[j] = 1;
  }
  static long nfail(0);
  for (size_t j(0); j < n; ++j)
    for (size_t l(j + 1); l < n; ++l) {
      if (!ok[j] || !ok[l] || !Beta2PairEligible((int)j, (int)l)) continue;
      const PhotonSet S{(int)j, (int)l};
      Vec4D_Vector legs, kp;
      Amplitude M2;
      PropShifts none;
      bool good(SubsetScaledLegs(S, legs, kp));
      Complex sj(0., 0.), sl(0., 0.);
      if (good) {
        Vec4D_Vector pp(legs);
        pp.push_back(kp[0]);
        pp.push_back(kp[1]);
        good = ComixRealNHel(pp, std::vector<int>{m_PhoHel[j], m_PhoHel[l]}, M2, none, true);
        if (good) {
          sj = StageEikonal(m_stage[j], legs, kp[0], m_PhoHel[j], (int)j);
          sl = StageEikonal(m_stage[l], legs, kp[1], m_PhoHel[l], (int)l);
          good = std::abs(sj) > 0. && std::abs(sl) > 0.;
        }
      }
      if (!good) {
        ++m_b2fail;
        ++nfail;
        long d(1);
        while (d*10 <= nfail) d *= 10;
        if (nfail == d)
          msg_Error()<<METHOD<<"(): no two-photon point or amplitude for "<<nfail
                     <<" pair(s) so far; their beta_2 is dropped."<<std::endl;
        continue;
      }
      ++m_bnamps[2];
      ++m_bnterms[2];
      Complex rest(1., 0.);
      for (size_t m(0); m < n; ++m) if (m != j && m != l) rest *= c[m];
      const Complex w(m_Sprod*rest), ss(sj*sl);
      const int tnpclass(m_tnp ? TnpSetClass(S) : tnp_b2m);
      for (int f(0); f < nh; ++f) {
        const Complex ujl(M2.m_A[f]/ss);
        const Complex add(w*(ujl - u[j].m_A[f]*c[l] - u[l].m_A[f]*c[j]
                             + c[j]*c[l]*B0.m_A[f]));
        if (IsBad(std::abs(add))) continue;
        m_AmpBeta2.m_A[f] += add;
        TnpAdd(tnpclass, f, add);
        if (m_treeprobe && m_treeorder.size() > 2) m_treeorder[2].m_A[f] += add;
      }
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
  // one radiating side (MODE ISR or FSR): M_n without the other side's
  // emission, as the generator's M_m (AddTreeBetas)
  PropShifts shifts;
  if (SingleSided()) AddPlacementVetoes(PhotonRange(n), shifts);
  const bool haveM(RealProvider(n) != NULL
                   && ComixRealNShifted(pp, PhotonRange(n), Mex, shifts, shifts.empty()));
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
