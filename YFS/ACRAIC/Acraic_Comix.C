#include "YFS/ACRAIC/Acraic_Base.H"
#include "YFS/NLO/NLO_Base.H"   // MapMomenta, for the reduced beta_1 kinematics
#include "YFS/NLO/Real_Correction.H"

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"

#include <cstdlib>
#include <iostream>

using namespace YFS;
using namespace ATOOLS;

namespace {

  /*!
    Flat Spin_Amplitudes index: bit i = leg i (0 -> +1, 1 -> -1), photons
    above the fermion legs. Same packing as Amplitude::Idx, so relabelling a
    set of legs is one FlipLegs.
  */
  inline size_t FlatIndex(size_t ferm, int hg, int nlegs)
  { return JoinLegs(ferm, (size_t)hg, nlegs); }


}




bool Acraic_Base::CalibrateComixMap(ComixCalib &c)
{
  if (p_bornproc == NULL || m_bornmomenta.size() < 4) return false;
  // The Born configuration, in its own rest frame.
  Vec4D_Vector bp(m_bornmomenta);
  Poincare bcms(bp[0] + bp[1]);
  for (size_t i(0); i < bp.size(); ++i) bcms.Boost(bp[i]);

  Amplitude cx, hand;
  const double sp_save(m_sp);
  m_sp = (bp[2] + bp[3]).Abs2();
  MakeProp(false);
  // MakePropT too: Comix's Born carries the t-channel at bp, so the hand Born
  // must as well (no-op except for Bhabha).
  MakePropT(bp);
  BornAmplitude(bp, hand);
  double cxme2(0.);
  const bool ok(ComixBornAmplitude(bp, cx, &cxme2));
  m_sp = sp_save;
  MakeProp(false);
  if (m_pacraic.size() >= 4) MakePropT(m_pacraic);
  if (!ok) return false;

  const int nh(Amplitude::NHel());
  double sh(0.), sc(0.);
  for (int f = 0; f < nh; ++f) {
    sh += std::norm(m_e*m_e*hand.m_A[f]);
    sc += std::norm(cx.m_A[f]);
  }
  if (!(sh > 0.) || !(sc > 0.)) return false;

  const double rh(1./sqrt(sh)), rc(1./sqrt(sc));
  int bestmask(-1);
  double bestmet(-1.), nextmet(-1.);
  for (int m = 0; m < nh; ++m) {
    double met(0.);
    for (int f = 0; f < nh; ++f)
      met += sqr(std::abs(m_e*m_e*hand.m_A[f])*rh - std::abs(cx.m_A[FlipLegs(f, m)])*rc);
    met = sqrt(met);
    if (bestmask < 0 || met < bestmet) { nextmet = bestmet;
                                         bestmet = met; bestmask = m; }
    else if (nextmet < 0. || met < nextmet) nextmet = met;
  }

  // |C/H| must be the same in every live slot. The threshold is coarse on
  // purpose: Comix fills the helicity-flip entries at O(m_e/E), the hand Born
  // sets them to zero.
  double rmin(0.), rmax(0.);
  int nlive(0);
  for (int f = 0; f < nh; ++f) {
    const Complex H(m_e*m_e*hand.m_A[f]), C(cx.m_A[FlipLegs(f, bestmask)]);
    if (std::abs(H) > 1e-3*sqrt(sh) && std::abs(C) > 1e-3*sqrt(sc)) {
      const double a(std::abs(C/H));
      if (nlive++ == 0) rmin = rmax = a;
      else { rmin = Min(rmin, a); rmax = Max(rmax, a); }
    }
  }

  c.mask = bestmask; c.met = bestmet; c.next = nextmet;
  c.N = sqrt(sc/sh);  c.sh = sh; c.sc = sc; c.me2 = cxme2;
  c.rmin = rmin; c.rmax = rmax; c.nlive = nlive;
  return true;
}


/*!
  Adopt the Born calibration for COMIX_REAL_FLIP < 0 and/or
  COMIX_REAL_NORM <= 0. Once per run: the map is a convention, not per event.
*/
void Acraic_Base::DeriveComixMap()
{
  m_comixcalibrated = true;
  const bool needflip(m_comixflip < 0), neednorm(!(m_comixnorm > 0.));
  // a normalisation set by hand in the card is taken as given
  if (!neednorm) m_normexact = true;
  if (!needflip && !neednorm) return;

  ComixCalib c;
  if (!CalibrateComixMap(c)) {
    // Refuse rather than fall back to values fitted for one process.
    msg_Error()<<METHOD<<"(): Born calibration failed; the Comix real cannot"
               <<" be normalised. Set ACRAIC: COMIX_REAL_NORM and"
               <<" COMIX_REAL_FLIP by hand.\n";
    return;
  }

  const size_t nl(m_flavs.size());
  // The photon bit sits above the s_nlegs helicity slots, not m_flavs.size()
  // (they differ when a leg is a scalar, e.g. H l+ l-).
  if (needflip) m_comixflip = c.mask | (m_comixphoflip ? Amplitude::NHel() : 0);

  // N = sqrt(initial spin states) exactly: the hand-coded amplitude is
  // spin-averaged (KKMC convention), Comix's sum_hel is not. The measured N
  // only has to agree.
  double nspin(1.);
  for (size_t i(0); i < 2 && i < nl; ++i) nspin *= m_flavs[i].IntSpin() + 1.;
  // YFS: USE_MODEL_ALPHA 0: hand Born at alpha(0), Comix at alpha_model, so
  // m_rescale_alpha = alpha(0)/alpha_model enters (1 otherwise).
  const double expected(nspin > 0. ? m_rescale_alpha/sqrt(nspin) : 0.);
  if (neednorm) {
    const bool agrees(expected > 0. && c.N > 0.
                      && std::abs(c.N*expected - 1.) < 1e-2);
    m_comixnorm = agrees ? expected : (c.N > 0. ? 1./c.N : 0.);
    m_normexact = agrees;
    if (!agrees)
      msg_Error()<<METHOD<<"(): measured N = "<<c.N<<" is not"
                 <<" sqrt(initial spin states) = "<<(expected>0.?1./expected:0.)
                 <<"; falling back to the measured value, but the amplitudes"
                 <<" are then not in the convention assumed here.\n";
  }

  msg_Debugging()<<METHOD<<"(): Comix -> ACRAIC map derived from the Born:\n"
            <<"  flip mask   = "<<m_comixflip<<"  (fermion bits "<<c.mask
            <<", metric "<<c.met<<" against "<<c.next<<" for the runner-up"
            <<(m_comixphoflip ? "; photon bit set" : "; photon bit clear")<<")\n"
            <<"  normalisation = "<<m_comixnorm<<"  (N = "<<c.N
            <<"; 1/sqrt(initial spin states) = "<<expected<<")\n"
            <<"  |C/H| over the "<<c.nlive<<" live helicities: ["
            <<c.rmin<<", "<<c.rmax<<"]\n";
}



/*!
  Comix's one-photon amplitude at a momentum set taken verbatim, mapped into
  ACRAIC's convention. The caller may hold the fermions fixed and move only
  the photon; such a set needs COMIX: MOMENTUM_PROJECTION: 0.
*/
bool Acraic_Base::ComixRealAt(const Vec4D_Vector &pp, int hel, Amplitude &M1,
                            const double propscale)
{
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL) return false;
  // Helicity legs (s_nlegs), not flavours: they differ for scalars (H l+ l-).
  const int nlg(Amplitude::s_nlegs);
  if (pp.size() != m_flavs.size() + 1) return false;
  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (prov->ComixAmplitudes(pp, propscale));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != NumLegSubsets<int>(nlg + 1)) return false;
  const int nh(Amplitude::NHel());
  const int hgi(FlipLegs(HelicityBit<int>(hel), ComixPhotonFlip(nlg)));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)), hgi, nlg)];
  return true;
}

bool Acraic_Base::ComixRealShifted(const Vec4D_Vector &pp, int hel,
                                 Amplitude &M1, const PropShifts &shifts)
{
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL) return false;
  const int nlg(Amplitude::s_nlegs);
  if (pp.size() != m_flavs.size() + 1) return false;
  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (prov->ComixAmplitudesShifts(pp, shifts));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != NumLegSubsets<int>(nlg + 1)) return false;
  const int nh(Amplitude::NHel());
  const int hgi(FlipLegs(HelicityBit<int>(hel), ComixPhotonFlip(nlg)));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)), hgi, nlg)];
  return true;
}

Acraic_Base::PropShifts Acraic_Base::StageShifts(int iphot) const
{
  PropShifts sh;
  if (m_stage.size() != m_allphotons.size()) return sh;
  // One entry per stage: (Comix mask, its photons in this partition). A W
  // decay stage's mask is the W's daughters (StageShiftMask).
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    const size_t mask(StageShiftMask((int)g));
    if (mask == 0) continue;
    sh.push_back(std::make_pair(mask, StagePhotonSum((int)g, iphot)));
  }
  AddResonanceEntries(sh);   // ACRAIC: RESONANCE_EXPONENTIATION (else none)
  return sh;
}

/*
  ACRAIC: TCHANNEL_SHIFT. Automatic: on for 2 -> 2 only, since LegsAt's
  "radiating pair" beyond 2 -> 2 is just the first two final fermions and
  the shift then drags single-W t-channel photons towards t = 0.
*/
bool Acraic_Base::ExchangeLineShiftsOn() const
{
  static const tristate::code mode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                   ["TCHANNEL_SHIFT"].Get<tristate::code>());
  if (mode == tristate::automatic) return m_flavs.size() == 4;
  return mode == tristate::on;
}

/*
  The partition's reduced legs as leg entries theta_i (p~_i - p_i) (theta = -1
  incoming, +1 outgoing), from which COMIX::Amplitude::SetPropShifts places
  every space-like exchange line at the partition's reduced invariant; only
  the pole moves, the spinors stay physical. iphot < 0: the Born
  (BornLegsAt(m_PXvec)); iphot >= 0: M_1 of that photon (PartitionLegs, all
  zero at n = 1).
*/
void Acraic_Base::AddExchangeLineShifts(int iphot, PropShifts &sh,
                                       const Vec4D_Vector *base) const
{
  if (!ExchangeLineShiftsOn() || m_pacraic.size() < 4
      || m_flavs.size() != m_pacraic.size()) return;
  Vec4D_Vector pb;
  const bool ok(iphot < 0 ? BornLegsAt(m_PXvec, pb) : PartitionLegs(iphot, pb));
  if (!ok || pb.size() != m_pacraic.size()) return;
  // base: the legs the amplitude is evaluated on (default the physical ones)
  const Vec4D_Vector &b(base && base->size() == m_pacraic.size() ? *base : m_pacraic);
  for (size_t i(0); i < m_pacraic.size(); ++i) {
    const Vec4D d((i < 2 ? -1. : 1.)*(pb[i] - b[i]));
    sh.push_back(std::make_pair(LegUnion(LegBit(i),
                                PHASIC::Process_Base::s_propshiftleg), d));
  }
}

bool Acraic_Base::Beta1DressedFsr() const
{
  static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BETA1_DRESSED_FSR"].Get<bool>());
  return on;
}

/*
  ACRAIC: BETA1_DRESSED_FSR. The two legs of final stage `stage` as they were
  before the stage's other photons (all but iphot; iphot < 0: all) were
  emitted: the pair is moved from Q to Q + K, the lepton directions kept in
  the pair's rest frame, the momenta rescaled to (Q + K)^2. Beams and every
  other leg unchanged. Two-leg final stages only.
*/
bool Acraic_Base::DressedStageLegs(int stage, int iphot, Vec4D_Vector &legs) const
{
  if (stage < 0 || stage >= (int)m_stagelegs.size() || m_stagereduces[stage]) return false;
  const std::vector<int> sl(StageExternalLegs(stage));
  if (sl.size() != 2 || sl[0] < 2 || sl[1] < 2 || sl[0] >= (int)m_pacraic.size()
      || sl[1] >= (int)m_pacraic.size()) return false;
  const Vec4D pa(m_pacraic[sl[0]]), pb(m_pacraic[sl[1]]);
  const Vec4D Q(pa + pb), Qd(Q + StagePhotonSum(stage, iphot));
  const double Q2(Q.Abs2()), Qd2(Qd.Abs2());
  const double ma(m_flavs[sl[0]].Mass()), mb(m_flavs[sl[1]].Mass());
  if (!(Q2 > 0.) || !(Q[0] > 0.) || !(Qd[0] > 0.) || !(Qd2 > sqr(ma + mb))) return false;
  Poincare toQ(Q), toQd(Qd);
  Vec4D a(pa);
  toQ.Boost(a);
  const Vec3D va(a);
  const double na(va.Abs());
  if (!(na > 0.)) return false;
  const Vec3D n(va/na);
  const double lam(0.5*sqrt(Max(0., sqr(Qd2 - ma*ma - mb*mb) - 4.*ma*ma*mb*mb)/Qd2));
  Vec4D qa(sqrt(lam*lam + ma*ma), lam*n), qb(sqrt(lam*lam + mb*mb), -lam*n);
  toQd.BoostBack(qa);
  toQd.BoostBack(qb);
  legs = m_pacraic;
  legs[sl[0]] = qa;
  legs[sl[1]] = qb;
  return true;
}

Acraic_Base::PropShifts Acraic_Base::DressedStageShifts(int iphot, int dressed,
                                                    const Vec4D &dshift,
                                                    const Vec4D_Vector &base) const
{
  PropShifts sh;
  if (m_stage.size() != m_allphotons.size()) return sh;
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    const size_t mask(StageShiftMask((int)g));
    if (mask == 0) continue;
    sh.push_back(std::make_pair(mask, (int)g == dressed ? dshift
                                                       : StagePhotonSum((int)g, iphot)));
  }
  AddResonanceEntries(sh);
  AddExchangeLineShifts(iphot, sh, &base);
  return sh;
}

bool Acraic_Base::BornHasExchangeLine()
{
  std::map<const PHASIC::Process_Base*, int>::const_iterator
    it(m_exchline.find(p_bornproc));
  if (it != m_exchline.end()) return it->second != 0;
  if (p_bornproc == NULL || m_pacraic.size() < 4
      || m_flavs.size() != m_pacraic.size()) return false;
  // Probe: without an exchange line a leg shift leaves the amplitude bit for
  // bit unchanged.
  Amplitude A0, A1;
  PropShifts none, leg;
  const double e(1e-3*m_pacraic[0][0]);
  for (size_t i(0); i < m_pacraic.size(); ++i)
    leg.push_back(std::make_pair(LegUnion(LegBit(i),
                                 PHASIC::Process_Base::s_propshiftleg),
                                 Vec4D(e, 0.3*e*(i+1), -0.2*e, 0.5*e)));
  if (!ComixBornShifted(m_pacraic, A0, none) || !ComixBornShifted(m_pacraic, A1, leg))
    return false;                          // not cached: try again next event
  int has(0);
  for (int f(0); f < Amplitude::NHel(); ++f)
    if (A0.m_A[f] != A1.m_A[f]) { has = 1; break; }
  m_exchline[p_bornproc] = has;
  msg_Info()<<"ACRAIC: Born "<<(p_bornproc ? p_bornproc->Name() : std::string("?"))
            <<(has ? " has" : " has no")<<" space-like exchange line"
            <<(has ? "s" : "")<<"."<<std::endl;
  return has != 0;
}

Vec4D Acraic_Base::PartitionShift(int iphot, bool reducing) const
{
  Vec4D d;
  if (m_stage.size() != m_allphotons.size()) return d;
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot
        && (m_stagereduces[m_stage[i]] != 0) == reducing) d += m_allphotons[i];
  return d;
}


/*
  A real phase-space point with beams at total X and the radiating pair at
  total Y, each back to back on shell; directions from the event, spectators
  untouched. Comix needs a real point, not just a moved propagator pole.
*/
bool Acraic_Base::LegsAt(const Vec4D &X, const Vec4D &Y, Vec4D_Vector &pb) const
{
  pb = m_pacraic;
  const double X2(X.Abs2());
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  if (X2 <= sqr(ma + mb)) return false;
  Poincare cmsX(X);
  // ACRAIC: REDUCED_AXIS: beams along z in X's frame, as the generator
  // builds its s' point (matters for space-like exchange lines); false uses
  // the boosted-beam axis.
  static const bool reducedaxis(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                ["REDUCED_AXIS"].SetDefault(true).Get<bool>());
  Vec3D da(0., 0., 1.);
  if (!reducedaxis) {
    Vec4D ra(m_pacraic[0]);
    cmsX.Boost(ra);
    const double rap(Vec3D(ra).Abs());
    if (!(rap > 0.)) return false;
    da = Vec3D(ra)/rap;
  }
  const double kala(sqr(X2 - ma*ma - mb*mb) - 4.*ma*ma*mb*mb);
  if (kala <= 0.) return false;
  const double pa(sqrt(kala)/(2.*sqrt(X2)));
  Vec4D b1(sqrt(pa*pa + ma*ma), pa*da), b2(sqrt(pa*pa + mb*mb), -pa*da);
  cmsX.BoostBack(b1); cmsX.BoostBack(b2);
  pb[0] = b1; pb[1] = b2;
  // the radiating pair balances Y, along the direction it has in Y's frame
  const double Y2(Y.Abs2());
  const double m3(m_flavs[m_if1].Mass()), m4(m_flavs[m_if2].Mass());
  if (Y2 <= sqr(m3 + m4)) return false;
  Poincare cmsY(Y);
  Vec4D r3(m_pacraic[m_if1]);
  cmsY.Boost(r3);
  const double r3p(Vec3D(r3).Abs());
  if (!(r3p > 0.)) return false;
  const Vec3D d3(Vec3D(r3)/r3p);
  const double kalf(sqr(Y2 - m3*m3 - m4*m4) - 4.*m3*m3*m4*m4);
  if (kalf <= 0.) return false;
  const double pf(sqrt(kalf)/(2.*sqrt(Y2)));
  Vec4D q3(sqrt(pf*pf + m3*m3), pf*d3), q4(sqrt(pf*pf + m4*m4), -pf*d3);
  cmsY.BoostBack(q3); cmsY.BoostBack(q4);
  pb[m_if1] = q3; pb[m_if2] = q4;
  return true;
}

bool Acraic_Base::GeneratorBornAt(const Vec4D &R, Vec4D_Vector &pb) const
{
  const size_t nl(m_flavs.size());
  if (nl < 4 || m_prefsr.size() != nl) return false;
  Vec4D Q;
  for (size_t i(2); i < nl; ++i) Q += m_prefsr[i];
  const double s2(R.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (!(s2 > sqr(m1 + m2)) || !(Q.Abs2() > 0.) || !(Q[0] > 0.)) return false;
  Poincare toQ(Q);
  std::vector<Vec3D> q; std::vector<double> mm2; double msum(0.);
  for (size_t i(2); i < nl; ++i) {
    Vec4D qi(m_prefsr[i]); toQ.Boost(qi); q.push_back(Vec3D(qi));
    mm2.push_back(sqr(m_flavs[i].Mass())); msum += m_flavs[i].Mass();
  }
  const double M(sqrt(s2));
  if (msum >= M) return false;
  auto etot = [&](double xi) { double e(0.);
    for (size_t j(0); j < q.size(); ++j) e += sqrt(mm2[j] + xi*xi*q[j].Sqr());
    return e; };
  double lo(0.), hi(1.);
  while (etot(hi) < M && hi < 1e6) hi *= 2.;
  for (int it(0); it < 200; ++it) { const double mid(0.5*(lo+hi)); (etot(mid) < M ? lo : hi) = mid; }
  const double xi(0.5*(lo+hi));
  pb.assign(nl, Vec4D());
  const double sgn(m_bornmomenta.size() > 0 && m_bornmomenta[0][3] < 0. ? -1. : 1.);
  const double lam(0.5*sqrt(Max(0., sqr(s2 - m1*m1 - m2*m2) - 4.*m1*m1*m2*m2)/s2));
  pb[0] = Vec4D(sqrt(lam*lam + m1*m1), 0., 0.,  sgn*lam);
  pb[1] = Vec4D(sqrt(lam*lam + m2*m2), 0., 0., -sgn*lam);
  for (size_t i(2); i < nl; ++i) {
    const Vec3D v(xi*q[i-2]);
    pb[i] = Vec4D(sqrt(mm2[i-2] + v.Sqr()), v);
  }
  return true;
}

// ACRAIC: CRUDE_BORN_LEGS stage: each final pair taken back to before its photons
bool Acraic_Base::StageBornAt(const Vec4D &R, Vec4D_Vector &pb) const
{
  const size_t nl(m_flavs.size());
  if (nl < 4 || m_pacraic.size() != nl || m_stage.size() != m_allphotons.size())
    return false;
  const double s2(R.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (!(s2 > sqr(m1 + m2)) || !(R[0] > 0.)) return false;
  Vec4D_Vector fin(m_pacraic.begin() + 2, m_pacraic.end());
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    if ((int)g == m_initstage) continue;
    const std::vector<int> legs(StageExternalLegs((int)g));
    const Vec4D K(StagePhotonSum((int)g, -1));
    if (legs.size() != 2 || legs[0] < 2 || legs[1] < 2 || legs[0] >= (int)nl
        || legs[1] >= (int)nl || !(K[0] > 0.)) continue;
    const Vec4D pa(m_pacraic[legs[0]]), pbm(m_pacraic[legs[1]]);
    const Vec4D Q(pa + pbm), Qd(Q + K);
    const double Q2(Q.Abs2()), Qd2(Qd.Abs2());
    const double ma(m_flavs[legs[0]].Mass()), mb(m_flavs[legs[1]].Mass());
    if (!(Q2 > 0.) || !(Q[0] > 0.) || !(Qd[0] > 0.) || !(Qd2 > sqr(ma + mb))) return false;
    Poincare toQ(Q), toQd(Qd);
    Vec4D a(pa);
    toQ.Boost(a);
    const Vec3D va(a);
    const double na(va.Abs());
    if (!(na > 0.)) return false;
    const Vec3D n(va/na);
    const double lam(0.5*sqrt(Max(0., sqr(Qd2 - ma*ma - mb*mb) - 4.*ma*ma*mb*mb)/Qd2));
    Vec4D qa(sqrt(lam*lam + ma*ma), lam*n), qb(sqrt(lam*lam + mb*mb), -lam*n);
    toQd.BoostBack(qa);
    toQd.BoostBack(qb);
    fin[legs[0] - 2] = qa;
    fin[legs[1] - 2] = qb;
  }
  // the dressed final state sums to R; into R's rest frame, beams along z
  Poincare toR(R);
  pb.assign(nl, Vec4D());
  const double sgn(m_bornmomenta.size() > 0 && m_bornmomenta[0][3] < 0. ? -1. : 1.);
  const double lam(0.5*sqrt(Max(0., sqr(s2 - m1*m1 - m2*m2) - 4.*m1*m1*m2*m2)/s2));
  pb[0] = Vec4D(sqrt(lam*lam + m1*m1), 0., 0.,  sgn*lam);
  pb[1] = Vec4D(sqrt(lam*lam + m2*m2), 0., 0., -sgn*lam);
  for (size_t i(2); i < nl; ++i) { Vec4D q(fin[i - 2]); toR.Boost(q); pb[i] = q; }
  return true;
}

bool Acraic_Base::BornLegsAt(const Vec4D &X, Vec4D_Vector &pb) const
{
  Vec4D Y(X);
  for (size_t i(2); i < m_pacraic.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pacraic[i];
  return LegsAt(X, Y, pb);
}

/*
  The legs for photon iphot in the current partition: the other photons the
  partition counts as initial-state come out of the beams (X_wp), the rest
  go into the radiating pair (X_wp - k). Physical point at one photon.
*/
bool Acraic_Base::PartitionLegs(int iphot, Vec4D_Vector &pb) const
{
  if (iphot < 0 || iphot >= (int)m_allphotons.size()
      || m_stage.size() != m_allphotons.size()) return false;
  Vec4D X(m_pacraic[0] + m_pacraic[1]);
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot && m_stagereduces[m_stage[i]]) X -= m_allphotons[i];
  Vec4D Y(X - m_allphotons[iphot]);
  for (size_t i(2); i < m_pacraic.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pacraic[i];
  return LegsAt(X, Y, pb);
}

tchmultiphoton::code Acraic_Base::TchannelMultiphotonMode() const
{
  static const tchmultiphoton::code mode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                         ["TCHANNEL_MULTIPHOTON"].Get<tchmultiphoton::code>());
  return mode;
}

/*
  ACRAIC: TCHANNEL_MULTIPHOTON >= 1. The one-photon point for photon iphot,
  as YFS.NLO's REAL_MAP 2 (NLO_Base::MapMomentaScaled): the event with only
  this photon, scaled by x = sqrt(Q^2/(P - k_j)^2) so the final state keeps
  Q^2, built in the P rest frame and carried back by the generator's boost
  Poincare(Q), so that its soft limit is exactly BornLegsAt(Q). Used instead
  of PartitionLegs, whose beams are tilted by a hard companion's pT.
  legs gets the Born legs, kpt the photon.
*/
bool Acraic_Base::OnePhotonScaledLegs(int iphot, Vec4D_Vector &legs,
                                    Vec4D &kpt) const
{
  if (iphot < 0 || iphot >= (int)m_allphotons.size()
      || m_pacraic.size() != m_flavs.size() || m_pacraic.size() < 4) return false;
  const Vec4D P(m_pacraic[0] + m_pacraic[1]);
  Vec4D Q;
  for (size_t i(2); i < m_pacraic.size(); ++i) Q += m_pacraic[i];
  const Vec4D R(P - m_allphotons[iphot]);
  const double P2(P.Abs2()), Q2(Q.Abs2()), R2(R.Abs2());
  if (!(P2 > 0.) || !(Q2 > 0.) || !(R2 > 0.) || !(Q[0] > 0.)) return false;
  const double x(Min(1., sqrt(Q2/R2)));
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  const double sx(x*x*P2);
  const double kal(sqr(sx - ma*ma - mb*mb) - 4.*ma*ma*mb*mb);
  if (!(sx > sqr(ma + mb)) || !(kal > 0.)) return false;
  // in the P rest frame: beams back to back along z at x sqrt(P^2), the
  // photon scaled by x, the final state carried from Q to Q' = x(P - k)
  Poincare toP(P), toQ(Q);
  Vec4D k(m_allphotons[iphot]);
  toP.Boost(k);
  k *= x;
  // leg 0 along +z, as LegsAt's reduced axis (REDUCED_AXIS 1) puts it
  const double pz(sqrt(kal)/(2.*sqrt(sx)));
  const Vec4D pa(sqrt(pz*pz + ma*ma), 0., 0.,  pz);
  const Vec4D pb(sqrt(pz*pz + mb*mb), 0., 0., -pz);
  const Vec4D Qp(pa + pb - k);
  if (!(Qp[0] > 0.) || !(Qp.Abs2() > 0.)) return false;
  Poincare toQp(Qp);
  legs = m_pacraic;
  legs[0] = pa; legs[1] = pb;
  for (size_t i(2); i < m_pacraic.size(); ++i) {
    Vec4D q(m_pacraic[i]);
    toQ.Boost(q);           // the generator's Born frame
    toQp.BoostBack(q);      // to Q' in the P rest frame
    legs[i] = q;
  }
  // back by the generator's boost: the soft limit is BornLegsAt(Q)
  for (Vec4D &l : legs) toQ.BoostBack(l);
  toQ.BoostBack(k);
  kpt = k;
  Vec4D bal(legs[0] + legs[1] - kpt);
  for (size_t i(2); i < legs.size(); ++i) bal -= legs[i];
  const double scale(Max(legs[0][0] + legs[1][0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-10*scale) return false;
  return true;
}

/*
  ACRAIC: TCHANNEL_MULTIPHOTON 2. Factorised form of the partition's beta_1
  sum a + sum_j v_j (a = m_Sprod beta_0), the amplitude-level analogue of
  YFS.NLO's REAL_COMBINE product, for Borns near a pole of the reduced frame:
      c_j = <a|a + v_j>/<a|a>,   r_j = v_j - (c_j - 1) a,
      T   = a prod_j c_j + sum_j r_j prod_{l != j} c_l .
  T = a + sum_j v_j to first order (exactly for one photon); the difference
  T - a - sum_j v_j is added.
*/
void Acraic_Base::AddFactorisedRemainder()
{
  const size_t nt(m_b1terms.size());
  if (nt < 2) return;
  const int nh(Amplitude::NHel());
  std::vector<Complex> a(nh);
  double na(0.);
  for (int f(0); f < nh; ++f) {
    a[f] = m_Sprod * m_partborn0.m_A[f];
    na += std::norm(a[f]);
  }
  if (!(na > 0.) || IsBad(na)) return;
  std::vector<Complex> c(nt, Complex(1., 0.));
  for (size_t j(0); j < nt; ++j) {
    Complex proj(0., 0.);
    for (int f(0); f < nh; ++f) proj += std::conj(a[f]) * m_b1terms[j].m_A[f];
    c[j] += proj/na;
  }
  auto product_without = [&](size_t skip) {
    Complex p(1., 0.);
    for (size_t l(0); l < nt; ++l) if (l != skip) p *= c[l];
    return p; };
  const Complex call(product_without(nt));
  for (int f(0); f < nh; ++f) {
    Complex additive(a[f]), factorised(a[f] * call);
    for (size_t j(0); j < nt; ++j) {
      const Complex vj(m_b1terms[j].m_A[f]);
      additive += vj;
      factorised += (vj - (c[j] - 1.) * a[f]) * product_without(j);
    }
    const Complex d(factorised - additive);
    if (IsBad(std::abs(d))) continue;
    m_AmpExpo1.m_A[f]    += d;
    m_AmpBornReal.m_A[f] += d;
    TnpAdd(tnp_multi, f, d);
    // ACRAIC: TNP_IFI: the factorised remainder taken as charge even
    TnpIfiAdd(f, Complex(0., 0.), m_tnppartflip*d);
  }
}

/*
  Total eikonal sum_l Q_l theta_l SfactorLeg(p_l) (eq. 6.4); equals
  sum_g m_Sfac[g][j] on the physical legs, recomputed here on reduced legs.
*/
Complex Acraic_Base::TotalEikonal(const Vec4D_Vector &p, const Vec4D &k,
                                int hel)
{
  Complex tot(0., 0.);
  for (size_t l(0); l < m_flavs.size() && l < p.size(); ++l) {
    // ISR only: the final current does not radiate (see m_isronly); FSR only:
    // the initial current does not (m_fsronly)
    if (m_isronly && l >= 2) continue;
    if (m_fsronly && l < 2) continue;
    const double w(m_flavs[l].Charge() * (l < 2 ? -1. : +1.));
    if (w != 0.) tot += w * SfactorLeg(p[l], k, hel);
  }
  return tot;
}

Complex Acraic_Base::StageEikonal(int stage, const Vec4D_Vector &p,
                                const Vec4D &k, int hel, int iphot)
{
  // every leg of the stage, a resonance built from its daughters in p
  return StageCurrent(stage, iphot, k, hel, p);
}

/*
  beta_1 for one photon from Comix, partition independent:

      beta_1(k) = M_1(k) - s_tot(k) B ,

  s_tot the eikonal of every leg, B the Born on the same legs; no stage
  split, pseudo-flux or (1 - CKine) needed, since

      M_n  ~  sum_wp prod s^wp beta_0(X_wp)  +  sum_j prod_{i!=j} s_tot(k_i)
                                                * beta_1(k_j).

  Adding it on one designated stage of photon j visits each assignment of
  the others once: sum_{wp: wp_j = 0} prod_{i!=j} s^{wp_i} = prod_{i!=j} s_tot(k_i).
*/
bool Acraic_Base::ComixInfraredSubtracted_1_0(const Vec4D &k, int hel,
                                            int iphot)
{
  if (!m_cxbalignok || iphot < 0) return false;
  if (iphot >= (int)m_stage.size() || (int)m_Sfac.size() != m_nstages
      || iphot >= (int)m_Sfac[m_stage[iphot]].size()) return false;
  /*
    M_1 goes on one partition per assignment of the others (photon on the
    reducing stage, never "index 0", which is the non-radiating stage for
    nu nu); the beta_0 subtraction goes in every partition. A collapsed
    (fixed-stage) photon is never enumerated, so it gets M_1 in every
    partition against the total eikonal.
  */
  const bool fixed(m_fixedstage.size() == m_stage.size()
                   && m_fixedstage[iphot] != 0);
  // YFS: IFI false: partitions summed incoherently, each carries its own
  // share of M_1 (AddPlacementVetoes).
  const bool addm1(!m_ifi || fixed || m_stagereduces[m_stage[iphot]]);
  const double rn(RealNorm());
  if (!(rn > 0.)) return false;
  /*
    Per assignment of the others, the two partitions differing only in this
    photon's placement give (hep-ph/0006359)

      prod_{i!=j} s_i^{wp_i} [ M_1^I(k_j; X_{j in I}) + M_1^F(k_j; X_{j in F})
                              - s_I(k_j) beta_0(X_{j in I})
                              - s_F(k_j) beta_0(X_{j in F}) ],

    M_1^{I,F} at the physical spinors with propagators at the partition's X.
    ACRAIC: BETA1_LEGS: 0 (default): physical legs + k_j with per-stage
    propagator shifts (ComixRealShifted). 1: balanced point from
    PartitionLegs (right propagators, wrong spinors), kept for comparison.
  */
  static const beta1legs::code legsmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                            ["BETA1_LEGS"].Get<beta1legs::code>());
  const bool rebuild(legsmode == beta1legs::balanced || m_redborn);
  Vec4D_Vector pb;
  Vec4D dI;
  PropShifts shifts;
  if (rebuild) { if (!PartitionLegs(iphot, pb)) return false; }
  else {
    pb = m_pacraic;
    shifts = StageShifts(iphot);
    AddExchangeLineShifts(iphot, shifts); // space-like exchange lines
    // YFS: IFI false, or one radiating side (ISR or FSR only)
    if (!m_ifi || SingleSided()) AddPlacementVetoes(std::vector<int>(1, iphot), shifts);
    dI = PartitionShift(iphot, true);   // the initial stage's photons, for S
  }
  // Softness y = 2k.S/S^2 against the radiating system S (beams less the
  // other initial-state photons); independent of this photon's placement, so
  // both partitions of one beta_1 term are cut together. Below ACRAIC:
  // BETA1_XCUT beta_1 = 0: no M_1, no subtraction.
  Vec4D S(rebuild ? pb[0] + pb[1] : m_pacraic[0] + m_pacraic[1] - dI);
  // A photon on a W decay stage: softness against P_W as produced
  // (daughters + the stage's photons).
  if (WStagesActive() && !rebuild && m_stage[iphot] != m_initstage
      && !fixed) {
    const int g(m_stage[iphot]);
    for (size_t l(0); l < m_stagelegs[g].size(); ++l)
      if (IsResonanceLeg(m_stagelegs[g][l].leg))
        S = StageSystemMomentum(g) + StagePhotonSum(g, -1);
  }
  const double S2(S.Abs2());
  const double y(S2 > 0. ? 2.*(k*S)/S2 : 0.);
  static const double xcut(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                           ["BETA1_XCUT"].Get<double>());
  if (!(y > xcut)) return true;     // m_b2on[iphot] stays 0: beta_1 = 0 here
  /*
    The subtraction is m_partborn0 (the Born InfraredSubtractedME_0_0 added,
    so the (1 - n) cancellation is exact) times the other photons' eikonals
    and this photon's stage eikonal on the legs Comix was handed.

    ACRAIC: TCHANNEL_MULTIPHOTON >= 1: M_1 at OnePhotonScaledLegs, used as
    R_j = M_1/s_j(point) times the physical m_sactu[j] (rescale), so the
    subtraction is exactly s_phys B_0.
  */
  Vec4D km(k);
  bool onept(false);
  if (addm1 && rebuild && TchannelMultiphotonActive()
      && m_sactu.size() == m_allphotons.size()) {
    Vec4D_Vector po;
    Vec4D ko;
    if (OnePhotonScaledLegs(iphot, po, ko)) { pb = po; km = ko; onept = true; }
    else {
      ++m_tmpfail;
      long d(1);
      while (d*10 <= m_tmpfail) d *= 10;
      if (m_tmpfail == d)
        msg_Error()<<METHOD<<"(): ACRAIC: TCHANNEL_MULTIPHOTON: no one-photon "
                   <<"point for "<<m_tmpfail<<" photon(s) so far; those keep "
                   <<"PartitionLegs."<<std::endl;
    }
  }
  /*
    ACRAIC: BETA1_DRESSED_FSR (2026-10-11). M_1(k_j) is computed once, in the
    partition with photon j on the initial stage, and carries the emission
    from the final pair too; the final pair it is handed is the physical one,
    which has already given its energy to the photons this partition puts on
    the final stage ("companions"). In deep radiative return Bhabha a lepton
    that lost 95% of its energy to collinear companions puts the other photons
    of the bunch inside its widened dead cone, where the helicity-flip part of
    M_1 is large: beta_1/(s beta_0) up to 1.5e3, event weights to 1e5.
    With companions, M_1 is taken at the final pair dressed with them
    (DressedStageLegs: the pair before they were emitted), with its own
    subtraction there, and the standard subtractions put back:
        M_1,eff = M_1(pd) - s_I(pd) B^I(pd) - s_F(pd) B^F(pd)
                  + s_I B_0^I + s_F B_0^F ,
    B^{I,F} the Born with photon j's momentum on the initial / final side
    lines (as m_partborn0 in the two partitions of j), s_{I,F} the stage
    eikonals, unadorned = physical legs. The code subtracts s_I B_0^I here and
    s_F B_0^F in the partner partition, so the net beta_1 is the dressed
    point's own, M_1(pd) - s_I(pd) B^I(pd) - s_F(pd) B^F(pd): zero as k_j
    goes soft, and the old one as the companions go soft (pd -> physical).
    Two stages (initial + one two-leg final) only; no companions: unchanged.
  */
  bool dressed(false);
  Amplitude M1fix;     // - s_I(pd) B^I(pd) - s_F(pd) B^F(pd) + s_F B_0^F (+ s_I B_0^I below)
  int gF(-1);
  if (addm1 && !rebuild && !fixed && m_ifi && Beta1DressedFsr() && !WStagesActive()
      && m_nstages == 2 && m_stage[iphot] == m_initstage
      && m_sactu.size() == m_allphotons.size()) {
    gF = 1 - m_initstage;
    bool companions(false);
    for (size_t i(0); i < m_allphotons.size(); ++i)
      if ((int)i != iphot && m_stage[i] == gF) { companions = true; break; }
    Vec4D_Vector pd;
    if (companions && DressedStageLegs(gF, iphot, pd)) {
      static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                              ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
      const double pf(pfmode == pseudoflux::rho0_and_rho1 ? m_pflux : 1.);
      const Vec4D rinit(StagePhotonSum(m_initstage, iphot));   // companions on I, not j
      const Vec4D rfin(StagePhotonSum(gF, iphot));             // the dressed-in ones
      // Born shifts: initial side by ri, final side by rf, space-like lines
      // at BornLegsAt(X) against the legs `base`
      auto bornshifts = [&](const Vec4D &ri, const Vec4D &rf, const Vec4D &X,
                            const Vec4D_Vector &base, PropShifts &sh) {
        sh.clear();
        for (size_t g(0); g < m_stagelegs.size(); ++g) {
          const size_t mask(StageShiftMask((int)g));
          if (mask == 0) continue;
          sh.push_back(std::make_pair(mask, (int)g == m_initstage ? ri : rf));
        }
        AddResonanceEntries(sh);
        Vec4D_Vector pbx;
        if (ExchangeLineShiftsOn() && m_pacraic.size() >= 4
            && m_flavs.size() == m_pacraic.size() && BornLegsAt(X, pbx)
            && pbx.size() == m_pacraic.size())
          for (size_t i(0); i < m_pacraic.size(); ++i)
            sh.push_back(std::make_pair(LegUnion(LegBit(i), PHASIC::Process_Base::s_propshiftleg),
                                        (i < 2 ? -1. : 1.)*(pbx[i] - base[i])));
      };
      const Vec4D XI(m_PXvec), XF(m_PXvec + k);   // j on I / on F
      PropShifts shI, shF, shF0;
      bornshifts(rinit + k, Vec4D(), XI, pd, shI);
      bornshifts(rinit, k, XF, pd, shF);
      bornshifts(rinit, rfin + k, XF, m_pacraic, shF0);
      Amplitude CI, CF, CF0;
      if (ComixBornShifted(pd, CI, shI) && ComixBornShifted(pd, CF, shF)
          && ComixBornShifted(m_pacraic, CF0, shF0)) {
        const Complex sIpd(StageEikonal(m_initstage, pd, k, hel, iphot));
        const Complex sFpd(StageEikonal(gF, pd, k, hel, iphot));
        const Complex sFph(StageEikonal(gF, m_pacraic, k, hel, iphot));
        const auto fl(ComixLegFlips(Amplitude::s_nlegs));
        for (int f(0); f < Amplitude::NHel(); ++f) {
          const Complex a(pf*m_cxbalign.m_A[f]);
          M1fix.m_A[f] = - sIpd*a*CI.m_A[FlipLegs(f, fl)] - sFpd*a*CF.m_A[FlipLegs(f, fl)]
                         + sFph*a*CF0.m_A[FlipLegs(f, fl)];
        }
        pb = pd;
        shifts = DressedStageShifts(iphot, gF, Vec4D(), pd);
        dressed = true;
      }
    }
    if (companions && !dressed) {
      ++m_dressfail;
      long d(1);
      while (d*10 <= m_dressfail) d *= 10;
      if (m_dressfail == d)
        msg_Error()<<METHOD<<"(): ACRAIC: BETA1_DRESSED_FSR: no dressed point for "
                   <<m_dressfail<<" photon(s) so far; those keep the physical legs."
                   <<std::endl;
    }
  }
  Amplitude M1;
  // ACRAIC: TNP_IFI: M_1^I, so that the flipped M_1 is M_1^I - M_1^F = 2 M_1^I - M_1
  Amplitude M1I;
  // with greco, M_1 and M_1^I of the flipped weight at a_g -> -a_g
  Amplitude M1b, M1Ib;
  bool tnpsplit(false);
  if (addm1) {
    Vec4D_Vector pp(pb);
    pp.push_back(km);
    const bool ok(rebuild ? ComixRealAt(pp, hel, M1, -1.)
                          : ComixRealShifted(pp, hel, M1, shifts));
    if (!ok) return false;
    if (m_tnpifi) {
      // the physical-leg, enumerated case only; elsewhere M_1 counts as even
      if (!rebuild && !fixed && !dressed && !onept && m_stage[iphot] == m_initstage) {
        tnpsplit = TnpIfiInitialM1(iphot, pp, hel, shifts, M1, M1I);
        M1b = M1; M1Ib = M1I;
        if (tnpsplit && !m_resexp.empty()) {
          const PropShifts shb(TnpIfiFlippedResonance(shifts));
          tnpsplit = ComixRealShifted(pp, hel, M1b, shb)
                     && TnpIfiInitialM1(iphot, pp, hel, shb, M1b, M1Ib);
        }
      }
      if (!tnpsplit) ++m_tnpifiunsplit;
    }
  }
  // ACRAIC: TNP_IFI: this photon's flip (1 on the production stage, where M_1 goes)
  const Complex tnpfj(m_tnpifi && iphot < (int)m_tnpflipph.size()
                      && std::abs(m_tnpflipph[iphot]) > 0. ? m_tnpflipph[iphot] : Complex(1., 0.));
  // the other photons' eikonal product; never sProd/s_j (diverges at zeros of s_j)
  Complex w(1., 0.);
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot)
      w *= (m_sactu.size() == m_allphotons.size() ? m_sactu[i] : m_Sfac[m_stage[i]][i]);
  Complex sj(fixed ? TotalEikonal(pb, km, hel)
                   : StageEikonal(m_stage[iphot], pb, km, hel, iphot));
  Complex rescale(1., 0.);
  if (onept) {
    if (!(std::abs(sj) > 0.)) return true;   // no eikonal at the point
    rescale = m_sactu[iphot]/sj;
    sj = m_sactu[iphot];
  }
  const int nh(Amplitude::NHel());
  // ACRAIC: ORDER 2: store this beta_1's M_1 and eikonal for beta_2's
  // subtraction (Acraic_Beta2.C).
  if (m_order == 2 && iphot < (int)m_b2on.size()) {
    m_b2on[iphot] = 1;
    if (iphot < (int)m_b2onept.size()) m_b2onept[iphot] = onept ? 1 : 0;
    m_b2sj[iphot] = sj;
    m_b2hasM1[iphot] = addm1 ? 1 : 0;
    if (addm1)
      for (int f(0); f < nh; ++f)
        m_b2M1[iphot].m_A[f] = rescale * m_cxbalign.m_A[f] * M1.m_A[f]/rn
          + (dressed ? M1fix.m_A[f] + sj*m_partborn0.m_A[f] : Complex(0., 0.));
  }
  // TCHANNEL_MULTIPHOTON 2: this photon's term, for AddFactorisedRemainder
  const bool keepterm(onept && TchannelMultiphotonMode() == tchmultiphoton::factorised);
  if (keepterm) m_b1terms.push_back(Amplitude());
  const int tnpclass(m_tnp ? TnpPhotonClass(iphot, fixed) : tnp_isr);
  for (int f(0); f < nh; ++f) {
    const Complex sub(w * sj * m_partborn0.m_A[f]);
    Complex v(-sub);
    if (addm1) {
      // BETA1_DRESSED_FSR: + the dressed point's subtractions and the standard
      // ones put back (sj = s_I on the physical legs here, B_0^I = m_partborn0)
      const Complex m1(rescale * w * m_cxbalign.m_A[f] * M1.m_A[f]/rn
                       + (dressed ? w*(M1fix.m_A[f] + sj*m_partborn0.m_A[f]) : Complex(0., 0.)));
      v += m1;
      if (iphot < (int)m_realphotM1.size()) m_realphotM1[iphot].m_A[f] += m1;
    }
    m_AmpExpo1.m_A[f]    += v;
    TnpAdd(tnpclass, f, v);
    if (m_tnpifi) {
      // the same term with the final-stage charges flipped
      Complex vb(-m_tnppartflip * w * sj * m_partborn0flip.m_A[f]);
      if (addm1) {
        const Complex m1b(tnpsplit ? rescale * w * m_cxbalign.m_A[f]
                                     * (2.*M1Ib.m_A[f] - M1b.m_A[f])/rn
                                   : rescale * w * m_cxbalign.m_A[f] * M1.m_A[f]/rn
                                     + (dressed ? w*(M1fix.m_A[f] + sj*m_partborn0.m_A[f])
                                                : Complex(0., 0.)));
        vb += m_tnppartflip/tnpfj*m1b;
      }
      TnpIfiAdd(f, Complex(0., 0.), vb);
    }
    m_AmpBornReal.m_A[f] += v;
    if (iphot < (int)m_realphot.size()) m_realphot[iphot].m_A[f] += v;
    m_snapReal.m_A[f]    += v;
    if (keepterm) m_b1terms.back().m_A[f] = v;
  }
  return true;
}
