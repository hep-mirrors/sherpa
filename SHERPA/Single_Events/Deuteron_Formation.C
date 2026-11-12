#include "SHERPA/Single_Events/Deuteron_Formation.H"
#include "HADRON_RESCATTERING/Main/Deuteron_Production.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace SHERPA;
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

namespace {
  const double pi = M_PI;
  // ALEPH antideuteron selection at the Z pole, arXiv:hep-ex/0604023
  const double aleph_pmin = 0.62, aleph_pmax = 1.03, aleph_cosmax = 0.95;
  // OPAL antideuteron search at the Z pole, momentum range of Eq. (21) in
  // arXiv:1504.07242
  const double opal_pmin = 0.35, opal_pmax = 1.1;
  // pair maps of Fig. 6 in arXiv:1504.07242: k vs |p1+p2|, both 0..5 GeV
  const size_t map_nbins = 50;
  const double map_max   = 5.;
  // ALEPH B2, arXiv:hep-ex/0604023 (Eqs. (4), (5)): antiprotons at half the
  // antideuteron momentum window, same angular range
  const double aleph_half_pmin = 0.5*aleph_pmin, aleph_half_pmax = 0.5*aleph_pmax;
  // ALICE multiplicity classes, arXiv:1902.09290 (Sec. 2.2, Table 1): V0M
  // estimator = charged particles in V0-A (2.8 <= eta <= 5.1) and V0-C
  // (-3.7 <= eta <= -1.7); INEL>0 = at least one charged particle in |eta| < 1
  const size_t mult_nmax = 600, mult_npt = 50;   // pT bins of 0.1 GeV up to 5 GeV
  const double mult_ptmax = 5.;
  const std::vector<std::string> mult_species = {"d","dbar","p","pbar"};
}

Deuteron_Formation::Deuteron_Formation() :
  m_done(false), m_nevents(0.), p_deuterons(nullptr)
{
  m_name = std::string("Deuteron_Formation");
  m_type = eph::Hadronization;
  // the settings block DEUTERON_FORMATION of the run card
  Scoped_Settings s{ Settings::GetMainSettings()["DEUTERON_FORMATION"] };
  m_on   = s["MODE"].SetDefault(false).Get<bool>();
  if (!m_on) return;
  std::string model = s["MODEL"].SetDefault("CrossSection").Get<std::string>();
  double inv_sigma0 = s["INV_SIGMA0"].SetDefault(1.8).Get<double>();   // 1/barn
  double b_coal     = s["B_COAL"].SetDefault(3.).Get<double>();        // fm
  // coalescence momentum, cut on k = |p_p - p_n| in the pair frame, and radius
  // for the impact parameter (B_COAL <= 0: momentum condition only)
  double p_coal     = s["P_COAL"].SetDefault(2.*0.197326/3.).Get<double>(); // GeV
  m_nsamp           = size_t(std::max(1,s["NSAMP"].SetDefault(10).Get<int>()));
  m_outdir          = s["OUTPUT"].SetDefault("DeuteronFormationFromXSecs").Get<std::string>();
  // Particles with lifetimes above 100 fm/c count as stable (paper, Sec. IV),
  // i.e. nucleons from decays of hadrons with Gamma < hbar c/(100 fm) are
  // excluded.
  m_widthcut = 0.1973269804/100.;
  CoalescenceModel cmodel = CoalescenceModel::CrossSectionModel;
  if      (model=="HardDisk")    cmodel = CoalescenceModel::HardDisk;
  else if (model=="Exponential") cmodel = CoalescenceModel::Exponential;
  else if (model!="CrossSection")
    THROW(fatal_error,"Unknown DEUTERON_FORMATION:MODEL '"+model+"'.");
  p_deuterons = new Deuteron_Coalescence(b_coal,p_coal,inv_sigma0,cmodel);
  m_name     += std::string(": ")+model;
  InitHistograms();
  if (cmodel!=CoalescenceModel::CrossSectionModel)
    m_name += ", p_coal = "+ToString(p_coal)+" GeV, b_coal = "+ToString(b_coal)+" fm";
  msg_Info()<<METHOD<<": model = "<<model<<", 1/sigma_0 = "<<inv_sigma0
            <<" 1/b, NSAMP = "<<m_nsamp<<", output to "<<m_outdir<<"/.\n";
  if (s["TEST"].SetDefault(false).Get<bool>()) p_deuterons->Tests(m_outdir);
}

Deuteron_Formation::~Deuteron_Formation() {
  for (Histogram_Map::iterator hit=m_histos.begin();hit!=m_histos.end();++hit)
    delete hit->second;
  if (p_deuterons) delete p_deuterons;
}

void Deuteron_Formation::InitHistograms() {
  // type 1: also accumulates sum of w^2, written as third column (error)
  for (const std::string & tag : {"formed","expected"}) {
    for (const std::string & sp : {"d","dbar"}) {
      std::string pre = tag+"_"+sp+"_";
      // 0.05 GeV bins, so that e.g. the BaBar bins are sums of bins
      m_histos[pre+"p"]       = new Histogram(1, 0.,5.,100);
      m_histos[pre+"pT"]      = new Histogram(1, 0.,5.,50);
      m_histos[pre+"y"]       = new Histogram(1,-5.,5.,50);
      // 1/(2 pi pT) d^2N/dpT dy for |y|<0.5, as measured by ALICE, in
      // 0.1 GeV bins as the ALICE data (arXiv:1709.08522)
      m_histos[pre+"invyield_midy"] = new Histogram(1,0.,5.,50);
    }
  }
  // (anti)protons entering the formation (no weak-decay products), as the
  // input spectra (cf. Fig. 8 of arXiv:1504.07242): invariant yield at |y|<0.5
  for (const std::string & sp : {"p","pbar"})
    m_histos["nucleon_"+sp+"_invyield_midy"] = new Histogram(1,0.,5.,50);
  // p + pbar at |y|<0.5 split by their immediate parent, see ProtonOrigin
  for (const std::string & origin : {"direct","Delta","Nstar","strangebaryon","other"})
    m_histos["proton_origin_"+origin+"_invyield_midy"] = new Histogram(1,0.,5.,50);
  // charged pions and kaons at |y|<0.5 (no weak-decay products), as the
  // ALICE pi/K/p spectra (arXiv:1504.00024): diagnostic only, to tell whether
  // a too hard proton spectrum is shared by the mesons
  for (const std::string & sp : {"pi_plus","pi_minus","K_plus","K_minus"})
    m_histos["meson_"+sp+"_invyield_midy"] = new Histogram(1,0.,10.,100);
  // diagnostics of the nucleons' partners, see FillPartners: for mid-rapidity
  // (anti)protons in bins of their pT, the smallest k to any (anti)neutron and
  // the impact parameter of that pair; the number of protons, of protons with
  // a neutron within k < 0.5 GeV, and the summed number of scatterings in
  // their events, all in proton pT; and the number of scatterings per event
  for (const std::string & bin : {"0.5-1","1-2","2-3","3-5"}) {
    m_histos["partner_knearest_pT"+bin] = new Histogram(1,0.,5.,100);
    m_histos["partner_bnearest_pT"+bin] = new Histogram(1,0.,100.,100);
  }
  m_histos["partner_protons_pT"]        = new Histogram(1,0.,5.,25);
  m_histos["partner_protons_k05_pT"]    = new Histogram(1,0.,5.,25);
  m_histos["partner_protons_nscatt_pT"] = new Histogram(1,0.,5.,25);
  // origin of the deuteron-forming pairs (same in PythiaReference/pythia_deuterons.cc), in
  // proton pT: number of neutrons within k < 0.5 GeV per proton; summed charged pT in a cone
  // DeltaR < 0.5 around the proton (jet activity), for all protons and for protons with a
  // partner; and whether proton and nearest partner are made directly or in hadron decays
  m_histos["partner_protons_npart05_pT"] = new Histogram(1,0.,5.,25);
  m_histos["partner_protons_cone_pT"]    = new Histogram(1,0.,5.,25);
  m_histos["partner_pairs_cone_pT"]      = new Histogram(1,0.,5.,25);
  for (const std::string & o : {"direct_direct","direct_decay","decay_decay"})
    m_histos["partner_pairs_"+o+"_pT"] = new Histogram(1,0.,5.,25);
  // antibaryons near the protons (same in PythiaReference/pythia_deuterons.cc), a test of
  // whether the pairs come from two nearby diquark pops (B ... Bbar ... B): in proton pT, the
  // protons and the protons with a partner that have a baryon of opposite baryon number within
  // k < 1 GeV; and for proton pT 1.5-5 GeV the k of the nearest such baryon
  m_histos["partner_protons_anti1_pT"] = new Histogram(1,0.,5.,25);
  m_histos["partner_pairs_anti1_pT"]   = new Histogram(1,0.,5.,25);
  m_histos["partner_protons_kanti_pT1.5-5"] = new Histogram(1,0.,5.,50);
  m_histos["partner_pairs_kanti_pT1.5-5"]   = new Histogram(1,0.,5.,50);
  m_histos["events_nscatt"]           = new Histogram(1,0.,60.,60);
  for (const std::string & pair : {"pn","pp","nn"})
    m_histos["pairs_k_"+pair] = new Histogram(1,0.,5.,100);
  for (const Deuteron_Channel & ch : Deuteron_Coalescence::Channels())
    m_histos["expected_k_"+Deuteron_Coalescence::Name(ch)] = new Histogram(1,0.,5.,100);
  // impact parameter of (anti)p n pairs [fm], used by HardDisk/Exponential
  m_histos["pairs_b_pn"] = new Histogram(1,0.,100.,100);
  for (const std::string & pair : {"pbarnbar","pbarpbar","nbarnbar","pn","pp","nn"})
    m_pairmaps[pair].assign(map_nbins*map_nbins,0.);
  m_mult_events.assign(mult_nmax,0.);
  m_mult_nch_sum.assign(mult_nmax,0.);
  for (const std::string & sp : mult_species) {
    m_mult_spectrum_sum[sp].assign(mult_nmax*mult_npt,0.);
    m_mult_spectrum_square_sum[sp].assign(mult_nmax*mult_npt,0.);
    m_mult_midy_count_sum[sp].assign(mult_nmax,0.);
    m_event_midy_spectrum[sp].assign(mult_npt,0.);
  }
}

// Diagnostic: the immediate parent of a (anti)proton.  "direct" (made in the
// hadronisation), "Delta" (Delta(1232) decay), "Nstar" (other non-strange
// baryon), "strangebaryon" (e.g. Sigma*, Lambda* strong decays), or "other"
// (mesons, heavy flavour).  Used to find out which protons differ between
// generators; it does not enter the deuteron formation.
std::string Deuteron_Formation::ProtonOrigin(Particle * part) const {
  Blob * blob = part->ProductionBlob();
  if (!blob || blob->Type()!=btp::Hadron_Decay || blob->NInP()==0) return "direct";
  const Flavour & mother = blob->InParticle(0)->Flav();
  if (!mother.IsBaryon()) return "other";
  long int kf = mother.Kfcode();
  if (kf==1114 || kf==2114 || kf==2214 || kf==2224) return "Delta";
  int q1 = (kf/1000)%10, q2 = (kf/100)%10, q3 = (kf/10)%10;
  if (q1>3 || q2>3 || q3>3) return "other";
  if (q1==3 || q2==3 || q3==3) return "strangebaryon";
  return "Nstar";
}

// True if the nucleon descends from a long-lived hadron (lifetime above
// 100 fm/c, e.g. a weak decay Lambda -> p pi-): such nucleons are produced far
// from the collision and are excluded from the formation (paper, Sec. IV).
// The decay chain is followed upwards through all hadron decays.
bool Deuteron_Formation::FromWeakDecay(Particle * part) const {
  Blob * blob = part->ProductionBlob();
  while (blob && blob->Type()==btp::Hadron_Decay && blob->NInP()>0) {
    Particle * mother = blob->InParticle(0);
    if (mother->Flav().Width()<m_widthcut) return true;
    blob = mother->ProductionBlob();
  }
  return false;
}

Return_Value::code Deuteron_Formation::Treat(Blob_List * bloblist) {
  // Nothing to do if the phase is switched off, this event has already been
  // treated, or the event is empty.
  if (!m_on || m_done || bloblist->empty()) return Return_Value::Nothing;
  if (!EventComplete(bloblist)) return Return_Value::Nothing;
  m_done = true;
  m_nevents += 1.;
  ResetEvent();

  size_t nV0 = 0, nch_mid = 0;
  bool   inel0 = false;
  CountEventActivity(bloblist,nV0,nch_mid,inel0);
  std::vector<Particle *> nucleons;
  CollectNucleons(bloblist,nucleons);
  FillMesons(bloblist);
  FillPartners(nucleons,bloblist);
  FillExpected(nucleons);
  bool formed = FormDeuterons(nucleons,bloblist);
  AddEventToRunSums(nV0,nch_mid,inel0);
  return formed ? Return_Value::Success : Return_Value::Nothing;
}

// The event is complete when no blob waits for showers, colour
// reconnections, hadronisation or hadron decays any more; otherwise a later
// pass of the event phases will come back here.  needs_hadronRescatter and
// needs_softUE are not checked: the first is never cleared when the
// rescattering phase is off, the second stays set on the beam blobs of
// minimum-bias events.
bool Deuteron_Formation::EventComplete(Blob_List * bloblist) const {
  // the blob status is a bit field: "pending" has one bit set per flag, and
  // Status() & pending is non-zero if the blob carries any of these flags
  const int pending = blob_status::needs_showers | blob_status::needs_harddecays |
    blob_status::needs_reconnections | blob_status::needs_hadronization |
    blob_status::needs_hadrondecays;
  for (Blob * blob : *bloblist) {
    if (blob->Status() & pending) return false;
  }
  return true;
}

void Deuteron_Formation::ResetEvent() {
  m_event_counts.clear();
  m_event_midy_count.clear();
  for (Vector_Map::iterator vit=m_event_midy_spectrum.begin();
       vit!=m_event_midy_spectrum.end();++vit)
    std::fill(vit->second.begin(),vit->second.end(),0.);
}

// Event activity for the ALICE multiplicity classes, from the final-state
// charged particles (Sherpa keeps particles with c tau > 1 cm stable, as
// ALICE's primary-particle definition): nV0 = number in the V0 acceptance,
// nch_mid = number at |eta|<0.5, inel0 = at least one at |eta|<1.
void Deuteron_Formation::CountEventActivity(Blob_List * bloblist,size_t & nV0,
                                            size_t & nch_mid,bool & inel0) const {
  for (Particle * part : bloblist->ExtractParticles(part_status::active,1)) {
    if (part->DecayBlob() || part->Flav().IntCharge()==0) continue;
    // four-momentum (E, px, py, pz) in the lab frame
    const Vec4D & mom = part->Momentum();
    if (Vec3D(mom).Abs()<=std::abs(mom[3])) continue;     // along the beam: no eta
    // pseudorapidity -ln tan(theta/2), theta = angle to the beam axis
    double eta = mom.Eta();
    if ((eta>=2.8 && eta<=5.1) || (eta>=-3.7 && eta<=-1.7)) nV0++;
    if (std::abs(eta)<0.5) nch_mid++;
    if (std::abs(eta)<1.)  inel0 = true;
  }
}

// All final-state (anti)nucleons that may form deuterons: protons and
// neutrons (and antiparticles) that have not decayed and do not stem from a
// weak decay.  The (anti)protons are also histogrammed.
void Deuteron_Formation::CollectNucleons(Blob_List * bloblist,
                                         std::vector<Particle *> & nucleons) {
  for (Particle * part : bloblist->ExtractParticles(part_status::active,1)) {
    if (!part->Flav().IsNucleon()) continue;          // only p, n, pbar, nbar
    // in-particles of e.g. rescattering blobs may still be flagged active
    if (part->DecayBlob()) { m_event_counts["skipped_nucleons_with_decayblob"] += 1.; continue; }
    if (FromWeakDecay(part)) continue;
    nucleons.push_back(part);
    if (part->Flav().Kfcode()==kf_p_plus) FillProton(part);
  }
}

// Final-state charged pions and kaons at |y|<0.5 not stemming from a weak
// decay, histogrammed as the protons.  Does not change the event.
// Diagnostic of the nucleon correlations that make deuterons.  For every
// (anti)proton at |y|<0.5: the (anti)neutron closest in k (the relative
// momentum of Deuteron_Coalescence::K), and the impact parameter of that pair;
// histogrammed in bins of the proton pT.  Also counts the parton-parton
// scatterings of the event (hard-collision blobs), a measure of how central
// the pp collision was.  Does not change the event.
void Deuteron_Formation::FillPartners(const std::vector<Particle *> & nucleons,
                                      Blob_List * bloblist) {
  size_t nscatt = 0;
  for (Blob * blob : *bloblist) {
    if (blob->Type()==btp::Hard_Collision) nscatt++;
  }
  m_histos["events_nscatt"]->Insert(double(nscatt),1.);
  // charged final-state particles for the cone sums around the protons
  std::vector<Particle *> charged;
  for (Particle * part : bloblist->ExtractParticles(part_status::active,1)) {
    if (part->DecayBlob() || part->Flav().IntCharge()==0) continue;
    const Vec4D & mom = part->Momentum();
    if (mom.PPerp()<0.15 || Vec3D(mom).Abs()<=std::abs(mom[3]) || std::abs(mom.Eta())>2.) continue;
    charged.push_back(part);
  }
  // baryons as the formation sees them: the long-lived ones (lifetime above 100 fm/c, e.g.
  // Lambda, decayed later by Sherpa) and the final-state ones, none made in a weak decay
  std::vector<Particle *> baryons;
  for (Blob * blob : *bloblist) {
    for (int i=0;i<blob->NOutP();++i) {
      Particle * part = blob->OutParticle(i);
      if (!part->Flav().IsHadron() || !part->Flav().IsBaryon() || FromWeakDecay(part)) continue;
      Blob * decay = part->DecayBlob();
      if (decay && (decay->Type()!=btp::Hadron_Decay || part->Flav().Width()>=m_widthcut)) continue;
      baryons.push_back(part);
    }
  }
  for (Particle * proton : nucleons) {
    if (proton->Flav().Kfcode()!=kf_p_plus) continue;
    const Vec4D & mom = proton->Momentum();
    double pT = mom.PPerp();
    if (std::abs(mom.Y())>=0.5 || pT<=0.) continue;
    // nearest (anti)neutron of the same baryon number
    Particle * nearest = NULL;
    double knearest = 1.e12;
    for (Particle * neutron : nucleons) {
      if (neutron->Flav().Kfcode()!=kf_n) continue;
      if (neutron->Flav().IsAnti()!=proton->Flav().IsAnti()) continue;
      double k = Deuteron_Coalescence::K(mom,neutron->Momentum());
      if (k<knearest) { knearest = k; nearest = neutron; }
    }
    m_histos["partner_protons_pT"]->Insert(pT,1.);
    m_histos["partner_protons_nscatt_pT"]->Insert(pT,double(nscatt));
    // number of partners within k < 0.5 GeV
    double npart = 0.;
    for (Particle * neutron : nucleons) {
      if (neutron->Flav().Kfcode()!=kf_n || neutron->Flav().IsAnti()!=proton->Flav().IsAnti()) continue;
      if (Deuteron_Coalescence::K(mom,neutron->Momentum())<0.5) npart += 1.;
    }
    m_histos["partner_protons_npart05_pT"]->Insert(pT,npart);
    // jet activity: summed charged pT within DeltaR < 0.5 of the proton (proton excluded)
    double cone = 0.;
    for (Particle * part : charged) {
      if (part==proton) continue;
      const Vec4D & q = part->Momentum();
      double dphi = std::abs(q.Phi()-mom.Phi());
      if (dphi>M_PI) dphi = 2.*M_PI-dphi;
      if (sqr(q.Eta()-mom.Eta())+sqr(dphi)<0.25) cone += q.PPerp();
    }
    m_histos["partner_protons_cone_pT"]->Insert(pT,cone);
    // nearest baryon of opposite baryon number
    double kanti = 1.e12;
    for (Particle * baryon : baryons) {
      if (baryon->Flav().IsAnti()==proton->Flav().IsAnti()) continue;
      kanti = std::min(kanti,Deuteron_Coalescence::K(mom,baryon->Momentum()));
    }
    bool highpT = pT>=1.5 && pT<5.;
    if (kanti<1.) m_histos["partner_protons_anti1_pT"]->Insert(pT,1.);
    if (highpT)   m_histos["partner_protons_kanti_pT1.5-5"]->Insert(kanti,1.);
    if (nearest==NULL) continue;
    if (knearest<0.5) {
      m_histos["partner_protons_k05_pT"]->Insert(pT,1.);
      m_histos["partner_pairs_cone_pT"]->Insert(pT,cone);
      if (kanti<1.) m_histos["partner_pairs_anti1_pT"]->Insert(pT,1.);
      if (highpT)   m_histos["partner_pairs_kanti_pT1.5-5"]->Insert(kanti,1.);
      bool pdirect = ProtonOrigin(proton)=="direct", ndirect = ProtonOrigin(nearest)=="direct";
      std::string o = (pdirect && ndirect) ? "direct_direct" : ((pdirect || ndirect) ? "direct_decay" : "decay_decay");
      m_histos["partner_pairs_"+o+"_pT"]->Insert(pT,1.);
    }
    std::string bin = "";
    if      (pT>=0.5 && pT<1.) bin = "0.5-1";
    else if (pT>=1.  && pT<2.) bin = "1-2";
    else if (pT>=2.  && pT<3.) bin = "2-3";
    else if (pT>=3.  && pT<5.) bin = "3-5";
    if (bin.empty()) continue;
    m_histos["partner_knearest_pT"+bin]->Insert(knearest,1.);
    double b = std::sqrt(Deuteron_Coalescence::ImpactParameter2(proton,nearest));
    m_histos["partner_bnearest_pT"+bin]->Insert(b,1.);
  }
}

void Deuteron_Formation::FillMesons(Blob_List * bloblist) {
  for (Particle * part : bloblist->ExtractParticles(part_status::active,1)) {
    if (part->DecayBlob()) continue;
    kf_code kf = part->Flav().Kfcode();
    if (kf!=kf_pi_plus && kf!=kf_K_plus) continue;
    const Vec4D & mom = part->Momentum();
    double pT = mom.PPerp(), y = mom.Y();
    if (std::abs(y)>=0.5 || pT<=0.) continue;
    if (FromWeakDecay(part)) continue;
    std::string sp = (kf==kf_pi_plus) ? "pi" : "K";
    sp += part->Flav().IsAnti() ? "_minus" : "_plus";
    m_histos["meson_"+sp+"_invyield_midy"]->Insert(pT,1./(2.*pi*pT));
  }
}

void Deuteron_Formation::FillProton(Particle * part) {
  const Vec4D & mom = part->Momentum();
  // transverse momentum sqrt(px^2+py^2) and rapidity 1/2 ln((E+pz)/(E-pz))
  double pT = mom.PPerp(), y = mom.Y();
  std::string sp = part->Flav().IsAnti() ? "pbar" : "p";
  // dy = 1 for |y|<0.5
  if (std::abs(y)<0.5 && pT>0.) {
    m_histos["nucleon_"+sp+"_invyield_midy"]->Insert(pT,1./(2.*pi*pT));
    // origin of the (anti)proton: its immediate parent
    std::string origin = ProtonOrigin(part);
    m_histos["proton_origin_"+origin+"_invyield_midy"]->Insert(pT,1./(2.*pi*pT));
    m_event_counts["proton_midy_from_"+origin] += 1.;
    Blob * prodblob = part->ProductionBlob();
    if (prodblob && prodblob->Type()==btp::Hadron_Decay && prodblob->NInP()>0) {
      // parents of antiprotons charge-conjugated, so that p and pbar add up
      Flavour parent = prodblob->InParticle(0)->Flav();
      if (part->Flav().IsAnti()) parent = parent.Bar();
      m_event_counts["proton_midy_parent_"+parent.IDName()] += 1.;
    }
    size_t ipt = size_t(pT/mult_ptmax*mult_npt);
    if (ipt<mult_npt) m_event_midy_spectrum[sp][ipt] += 1./(2.*pi*pT);
    m_event_midy_count[sp] += 1.;
  }
  m_event_counts["nucleon_"+sp+"_total"] += 1.;
  double p = Vec3D(mom).Abs();
  if (p>aleph_half_pmin && p<aleph_half_pmax && std::abs(mom[3])/p<aleph_cosmax)
    m_event_counts["nucleon_"+sp+"_ALEPH_halfwindow"] += 1.;
}

// Expected yields: every pair, every open channel, weight P_i, with the
// deuteron kinematics sampled NSAMP times.  Also fills the pair histograms.
void Deuteron_Formation::FillExpected(const std::vector<Particle *> & nucleons) {
  for (size_t i=0;i<nucleons.size();i++) {
    for (size_t j=i+1;j<nucleons.size();j++) {
      Particle * A = nucleons[i], * B = nucleons[j];
      int type = Deuteron_Coalescence::PairType(A,B);
      if (type==0) continue;                           // nucleon + antinucleon
      bool anti = type<0;
      Vec4D pA = A->Momentum(), pB = B->Momentum();
      Poincare boost(pA+pB);
      Vec4D pAcm = pA;
      boost.Boost(pAcm);
      // k = |pA - pB| in the pair c.m. frame, as in the paper
      double k = 2.*Vec3D(pAcm).Abs();
      size_t np = (A->Flav().Kfcode()==kf_p_plus)+(B->Flav().Kfcode()==kf_p_plus);
      std::string pair = np==1 ? "pn" : (np==2 ? "pp" : "nn");
      m_histos["pairs_k_"+pair]->Insert(k,1.);
      m_event_counts["pairs_"+pair] += 1.;
      std::string mapname = anti ? (np==1 ? "pbarnbar" : (np==2 ? "pbarpbar" : "nbarnbar")) : pair;
      size_t ip = size_t(Vec3D(pA+pB).Abs()/map_max*map_nbins);
      size_t ik = size_t(k/map_max*map_nbins);
      if (ip<map_nbins && ik<map_nbins) m_pairmaps[mapname][ik*map_nbins+ip] += 1.;
      double b2 = Deuteron_Coalescence::ImpactParameter2(A,B);
      if (np==1) m_histos["pairs_b_pn"]->Insert(std::sqrt(b2),1.);
      for (const Channel_Probability & cp : p_deuterons->Probabilities(A,B,b2)) {
        std::string name = Deuteron_Coalescence::Name(cp.channel);
        m_histos["expected_k_"+name]->Insert(k,cp.prob);
        m_event_counts[std::string("expected_")+(anti ? "dbar_" : "d_")+name] += cp.prob;
        for (size_t n=0;n<m_nsamp;n++) {
          std::vector<Vec4D> moms;
          if (!p_deuterons->FinalState(cp.channel,pA,pB,moms)) continue;
          FillDeuteron("expected",anti,moms[0],cp.prob/double(m_nsamp));
        }
      }
    }
  }
}

// Actual formation (paper, Sec. III F): each pair binds with probability P_i;
// a nucleon that has formed a deuteron is not used again.  The deuterons are
// added to the event.  Returns true if at least one deuteron was formed.
bool Deuteron_Formation::FormDeuterons(const std::vector<Particle *> & nucleons,
                                       Blob_List * bloblist) {
  bool formed = false;
  std::vector<bool> used(nucleons.size(),false);
  for (size_t i=0;i<nucleons.size();i++) {
    if (used[i]) continue;
    for (size_t j=i+1;j<nucleons.size();j++) {
      if (used[j]) continue;
      Particle * A = nucleons[i], * B = nucleons[j];
      if (Deuteron_Coalescence::PairType(A,B)==0) continue;
      Deuteron_Channel channel =
        p_deuterons->SelectChannel(A,B,Deuteron_Coalescence::ImpactParameter2(A,B));
      if (channel==Deuteron_Channel::none) continue;
      Blob * blob = p_deuterons->BuildBlob(channel,A,B);
      if (!blob) continue;
      bloblist->push_back(blob);
      used[i] = used[j] = formed = true;
      Particle * deut = blob->OutParticle(0);
      bool anti = deut->Flav().IsAnti();
      FillDeuteron("formed",anti,deut->Momentum(),1.);
      m_event_counts[std::string("formed_")+(anti ? "dbar_" : "d_")+
                     Deuteron_Coalescence::Name(channel)] += 1.;
      msg_Out()<<"   * "<<(anti ? "anti" : "")<<"deuteron formed via "
               <<Deuteron_Coalescence::Name(channel)<<", p = "<<deut->Momentum()<<"\n";
      break;                                  // nucleon i is used: next i
    }
  }
  return formed;
}

void Deuteron_Formation::FillDeuteron(const std::string & tag,const bool & anti,
                                      const Vec4D & mom,const double & weight) {
  std::string pre = tag+"_"+(anti ? "dbar" : "d")+"_";
  double p = Vec3D(mom).Abs(), pT = mom.PPerp(), y = mom.Y();
  m_histos[pre+"p"]->Insert(p,weight);
  m_histos[pre+"pT"]->Insert(pT,weight);
  m_histos[pre+"y"]->Insert(y,weight);
  // dy = 1 for |y|<0.5
  if (std::abs(y)<0.5 && pT>0.) {
    m_histos[pre+"invyield_midy"]->Insert(pT,weight/(2.*pi*pT));
    if (tag=="expected") {
      std::string sp = anti ? "dbar" : "d";
      size_t ipt = size_t(pT/mult_ptmax*mult_npt);
      if (ipt<mult_npt) m_event_midy_spectrum[sp][ipt] += weight/(2.*pi*pT);
      m_event_midy_count[sp] += weight;
    }
  }
  m_event_counts[pre+"total"] += weight;
  if (p>aleph_pmin && p<aleph_pmax && std::abs(mom[3])/p<aleph_cosmax)
    m_event_counts[pre+"ALEPH_window"] += weight;
  if (p>opal_pmin && p<opal_pmax) m_event_counts[pre+"OPAL_window"] += weight;
}

// Adds the numbers of this event to the run sums: every per-event count and
// its square (for the errors in summary.dat), and, for INEL>0 events, the
// multiplicity-class sums in the row of the event's V0 multiplicity.
void Deuteron_Formation::AddEventToRunSums(const size_t & nV0,const size_t & nch_mid,
                                           const bool & inel0) {
  for (Count_Map::const_iterator cit=m_event_counts.begin();
       cit!=m_event_counts.end();++cit) {
    m_count_sums[cit->first]        += cit->second;
    m_count_square_sums[cit->first] += cit->second*cit->second;
  }
  if (!inel0) return;
  size_t row = std::min(nV0,mult_nmax-1);
  m_mult_events[row]  += 1.;
  m_mult_nch_sum[row] += double(nch_mid);
  for (Vector_Map::const_iterator vit=m_event_midy_spectrum.begin();
       vit!=m_event_midy_spectrum.end();++vit) {
    const std::string & sp = vit->first;
    for (size_t i=0;i<mult_npt;i++) {
      m_mult_spectrum_sum[sp][row*mult_npt+i]        += vit->second[i];
      m_mult_spectrum_square_sum[sp][row*mult_npt+i] += vit->second[i]*vit->second[i];
    }
    m_mult_midy_count_sum[sp][row] += m_event_midy_count[sp];
  }
}

void Deuteron_Formation::CleanUp(const size_t & mode) { m_done = false; }

// Writes all results to the OUTPUT directory at the end of the run.
void Deuteron_Formation::Finish(const std::string &) {
  if (!m_on || m_nevents<=0.) return;
  MakeDir(m_outdir,true);
  WriteHistograms();
  WritePairMaps();
  WriteMultiplicityClasses();
  WriteSummary();
  msg_Info()<<METHOD<<": "<<m_nevents<<" events, results in "<<m_outdir<<"/.\n";
}

// <name>.dat: every histogram, divided by the number of events and the bin
// width, i.e. per event and per unit of the histogrammed variable.
void Deuteron_Formation::WriteHistograms() const {
  for (Histogram_Map::const_iterator hit=m_histos.begin();hit!=m_histos.end();++hit) {
    hit->second->Scale(1./(m_nevents*hit->second->BinSize()));
    hit->second->Output(m_outdir+"/"+hit->first+".dat");
  }
}

// pairmap_<pair>.dat: one line per (k, |p1+p2|) bin with the two bin
// centres and the number of pairs per event in that bin.
void Deuteron_Formation::WritePairMaps() const {
  double width = map_max/double(map_nbins);
  for (Vector_Map::const_iterator mit=m_pairmaps.begin();mit!=m_pairmaps.end();++mit) {
    My_Out_File out(m_outdir+"/pairmap_"+mit->first+".dat");
    out.Open();
    *out<<"# |p1+p2| [GeV]  k [GeV]  pairs per event per bin ("
        <<width<<" GeV x "<<width<<" GeV)\n";
    for (size_t ik=0;ik<map_nbins;ik++)
      for (size_t ip=0;ip<map_nbins;ip++)
        *out<<(ip+0.5)*width<<" "<<(ik+0.5)*width<<" "
            <<mit->second[ik*map_nbins+ip]/m_nevents<<"\n";
    out.Close();
  }
}

// Multiplicity classes.  One row per value of the V0 multiplicity nV0 that
// occurred (last row: overflow); the numbers are sums over the events of the
// row, not divided by the number of events, so that the analysis can combine
// rows into classes.
//   multiplicity_events.dat: nV0, events, sum of dNch(|eta|<0.5), and the
//                            summed dN/dy(|y|<0.5) of d, dbar, p, pbar
//   multiplicity_<sp>_invyield_midy.dat:    summed spectra, 50 pT bins
//   multiplicity_<sp>_invyield_midy_sq.dat: summed squares (for the errors)
void Deuteron_Formation::WriteMultiplicityClasses() const {
  My_Out_File out(m_outdir+"/multiplicity_events.dat");
  out.Open();
  *out<<"# INEL>0 events binned in the V0M charged multiplicity nV0 (V0-A 2.8<=eta<=5.1, "
      <<"V0-C -3.7<=eta<=-1.7; arXiv:1902.09290)\n"
      <<"# nV0  events  sum_dNch(|eta|<0.5)  sum_dN/dy(|y|<0.5): d  dbar  p  pbar\n";
  for (size_t row=0;row<mult_nmax;row++) {
    if (m_mult_events[row]<=0.) continue;
    *out<<row<<" "<<m_mult_events[row]<<" "<<m_mult_nch_sum[row];
    for (const std::string & sp : mult_species) *out<<" "<<m_mult_midy_count_sum.at(sp)[row];
    *out<<"\n";
  }
  out.Close();
  double binwidth = mult_ptmax/double(mult_npt);
  for (const std::string & sp : mult_species) {
    for (bool squares : {false,true}) {
      const std::vector<double> & vals =
        squares ? m_mult_spectrum_square_sum.at(sp) : m_mult_spectrum_sum.at(sp);
      My_Out_File o(m_outdir+"/multiplicity_"+sp+"_invyield_midy"+(squares ? "_sq" : "")+".dat");
      o.Open();
      *o<<"# rows: nV0 (as multiplicity_events.dat); columns: "<<mult_npt<<" pT bins of "<<binwidth
        <<" GeV from 0; "<<(squares ? "sum over events of the squared per-event values"
                                    : "sum over events of 1/(2 pi pT) dN/dpT dy, |y|<0.5")
        <<" [GeV^-2"<<(squares ? ", squared" : "")<<"]\n";
      // per unit pT (the per-event values were filled per bin)
      double norm = squares ? binwidth*binwidth : binwidth;
      for (size_t row=0;row<mult_nmax;row++) {
        if (m_mult_events[row]<=0.) continue;
        *o<<row;
        for (size_t i=0;i<mult_npt;i++) *o<<" "<<vals[row*mult_npt+i]/norm;
        *o<<"\n";
      }
      o.Close();
    }
  }
}

// summary.dat: for every per-event count, its mean over the run and the
// error of the mean, sqrt(variance/N), from the sums and sums of squares.
void Deuteron_Formation::WriteSummary() const {
  My_Out_File summary(m_outdir+"/summary.dat");
  summary.Open();
  *summary<<"# Deuteron_Formation: "<<m_name<<", 1/sigma_0 = "
          <<p_deuterons->InvSigma0()<<" 1/b, NSAMP = "<<m_nsamp
          <<", events = "<<m_nevents<<"\n"
          <<"# ALEPH_window: 0.62 < p < 1.03 GeV, |cos(theta)| < 0.95; "
          <<"ALEPH: dbar = (5.9 +- 1.8 +- 0.5)e-6 per hadronic Z event\n"
          <<"# key  per-event-mean  error\n";
  for (Count_Map::const_iterator cit=m_count_sums.begin();cit!=m_count_sums.end();++cit) {
    double mean = cit->second/m_nevents;
    double var  = std::max(0.,m_count_square_sums.at(cit->first)/m_nevents-mean*mean);
    *summary<<cit->first<<"  "<<mean<<"  "<<sqrt(var/m_nevents)<<"\n";
  }
  summary.Close();
}
