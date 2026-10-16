#include "SHERPA/Single_Events/Deuteron_Formation.H"
#include "HADRON_RESCATTERING/Main/Deuteron_Production.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/MyStrStream.H"
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
}

Deuteron_Formation::Deuteron_Formation() :
  m_done(false), m_nevents(0.), p_deuterons(nullptr)
{
  m_name = std::string("Deuteron_Formation");
  m_type = eph::Hadronization;
  auto s = Settings::GetMainSettings()["DEUTERON_FORMATION"];
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
  for (auto & h : m_histos) delete h.second;
  if (p_deuterons) delete p_deuterons;
}

void Deuteron_Formation::InitHistograms() {
  // type 1: also accumulates sum of w^2, written as third column (error)
  for (std::string tag : {"formed","expected"}) {
    for (std::string sp : {"d","dbar"}) {
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
  // input spectra (cf. Fig. 8 of arXiv:1504.07242): invariant yield at |gamma|<0.5
  for (std::string sp : {"p","pbar"})
    m_histos["nucleon_"+sp+"_invyield_midy"] = new Histogram(1,0.,5.,50);
  for (std::string pair : {"pn","pp","nn"})
    m_histos["pairs_k_"+pair] = new Histogram(1,0.,5.,100);
  for (const auto & ch : Deuteron_Coalescence::Channels())
    m_histos["expected_k_"+Deuteron_Coalescence::Name(ch)] = new Histogram(1,0.,5.,100);
  // impact parameter of (anti)p n pairs [fm], used by HardDisk/Exponential
  m_histos["pairs_b_pn"] = new Histogram(1,0.,100.,100);
  for (std::string pair : {"pbarnbar","pbarpbar","nbarnbar","pn","pp","nn"})
    m_pairmaps[pair].assign(map_nbins*map_nbins,0.);
  m_mult_events.assign(mult_nmax,0.);
  m_mult_nch.assign(mult_nmax,0.);
  for (std::string sp : {"d","dbar","p","pbar"}) {
    m_mult_yield[sp].assign(mult_nmax*mult_npt,0.);
    m_mult_yield2[sp].assign(mult_nmax*mult_npt,0.);
    m_mult_dndy[sp].assign(mult_nmax,0.);
    m_ev_yield[sp].assign(mult_npt,0.);
  }
}

bool Deuteron_Formation::FromWeakDecay(Particle * part) const {
  Blob * blob = part->ProductionBlob();
  while (blob && blob->Type()==btp::Hadron_Decay && blob->NInP()>0) {
    Particle * mother = blob->InParticle(0);
    if (mother->Flav().Width()<m_widthcut) return true;
    blob = mother->ProductionBlob();
  }
  return false;
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
      if (ipt<mult_npt) m_ev_yield[sp][ipt] += weight/(2.*pi*pT);
      m_ev_dndy[sp] += weight;
    }
  }
  m_evcounts[pre+"total"] += weight;
  if (p>aleph_pmin && p<aleph_pmax && std::abs(mom[3])/p<aleph_cosmax)
    m_evcounts[pre+"ALEPH_window"] += weight;
  if (p>opal_pmin && p<opal_pmax) m_evcounts[pre+"OPAL_window"] += weight;
}

Return_Value::code Deuteron_Formation::Treat(Blob_List * bloblist) {
  // Phases restart from the beginning after every Success, so this is
  // reached only once all hadron decays are done.  Run once per event.
  if (!m_on || m_done || bloblist->empty()) return Return_Value::Nothing;
  m_done = true;
  m_nevents += 1.;
  m_evcounts.clear();
  m_ev_dndy.clear();
  for (auto & ey : m_ev_yield) std::fill(ey.second.begin(),ey.second.end(),0.);

  // event activity for the ALICE multiplicity classes (final-state charged
  // particles; Sherpa keeps particles with c tau > 1 cm stable, as ALICE's
  // primary-particle definition)
  size_t nV0 = 0, nch_mid = 0;
  bool   inel0 = false;
  for (Particle * part : bloblist->ExtractParticles(part_status::active,1)) {
    if (part->DecayBlob() || part->Flav().IntCharge()==0) continue;
    const Vec4D & mom = part->Momentum();
    if (Vec3D(mom).Abs()<=std::abs(mom[3])) continue;     // no pT: no eta
    double eta = mom.Eta();
    if ((eta>=2.8 && eta<=5.1) || (eta>=-3.7 && eta<=-1.7)) nV0++;
    if (std::abs(eta)<0.5) nch_mid++;
    if (std::abs(eta)<1.)  inel0 = true;
  }

  std::vector<Particle *> nucleons;
  for (Particle * part : bloblist->ExtractParticles(part_status::active,1)) {
    if (!part->Flav().IsNucleon()) continue;
    // in-particles of e.g. rescattering blobs may still be flagged active
    if (part->DecayBlob()) { m_evcounts["skipped_nucleons_with_decayblob"] += 1.; continue; }
    if (FromWeakDecay(part)) continue;
    nucleons.push_back(part);
    if (part->Flav().Kfcode()==kf_p_plus) {
      const Vec4D & mom = part->Momentum();
      double pT = mom.PPerp(), y = mom.Y();
      std::string sp = part->Flav().IsAnti() ? "pbar" : "p";
      // dy = 1 for |y|<0.5
      if (std::abs(y)<0.5 && pT>0.) {
        m_histos["nucleon_"+sp+"_invyield_midy"]->Insert(pT,1./(2.*pi*pT));
        size_t ipt = size_t(pT/mult_ptmax*mult_npt);
        if (ipt<mult_npt) m_ev_yield[sp][ipt] += 1./(2.*pi*pT);
        m_ev_dndy[sp] += 1.;
      }
      m_evcounts["nucleon_"+sp+"_total"] += 1.;
      double p = Vec3D(mom).Abs();
      if (p>aleph_half_pmin && p<aleph_half_pmax && std::abs(mom[3])/p<aleph_cosmax)
        m_evcounts["nucleon_"+sp+"_ALEPH_halfwindow"] += 1.;
    }
  }

  // Expected yields: every pair, every open channel, weight P_i, with the
  // deuteron kinematics sampled NSAMP times.
  for (size_t i=0;i<nucleons.size();i++) {
    for (size_t j=i+1;j<nucleons.size();j++) {
      Particle * A = nucleons[i], * B = nucleons[j];
      int type = Deuteron_Coalescence::PairType(A,B);
      if (type==0) continue;
      bool anti = type<0;
      Vec4D pA = A->Momentum(), pB = B->Momentum();
      Poincare boost(pA+pB);
      Vec4D pAcm = pA;
      boost.Boost(pAcm);
      // k = |pA - pB| in the pair c.m. frame, as in the paper
      double k = 2.*Vec3D(pAcm).Abs();
      size_t np = (A->Flav().Kfcode()==kf_p_plus)+(B->Flav().Kfcode()==kf_p_plus);
      m_histos[np==1 ? "pairs_k_pn" : (np==2 ? "pairs_k_pp" : "pairs_k_nn")]->Insert(k,1.);
      m_evcounts[np==1 ? "pairs_pn" : (np==2 ? "pairs_pp" : "pairs_nn")] += 1.;
      std::string pname = np==1 ? "pn" : (np==2 ? "pp" : "nn");
      if (anti) pname = np==1 ? "pbarnbar" : (np==2 ? "pbarpbar" : "nbarnbar");
      size_t ip = size_t(Vec3D(pA+pB).Abs()/map_max*map_nbins);
      size_t ik = size_t(k/map_max*map_nbins);
      if (ip<map_nbins && ik<map_nbins) m_pairmaps[pname][ik*map_nbins+ip] += 1.;
      double b2 = Deuteron_Coalescence::ImpactParameter2(A,B);
      if (np==1) m_histos["pairs_b_pn"]->Insert(std::sqrt(b2),1.);
      for (const auto & cp : p_deuterons->Probabilities(A,B,b2)) {
        std::string name = Deuteron_Coalescence::Name(cp.channel);
        m_histos["expected_k_"+name]->Insert(k,cp.prob);
        m_evcounts[std::string("expected_")+(anti ? "dbar_" : "d_")+name] += cp.prob;
        for (size_t n=0;n<m_nsamp;n++) {
          std::vector<Vec4D> moms;
          if (!p_deuterons->FinalState(cp.channel,pA,pB,moms)) continue;
          FillDeuteron("expected",anti,moms[0],cp.prob/double(m_nsamp));
        }
      }
    }
  }

  // Actual formation; nucleons used once are excluded from further pairs.
  bool formed = false;
  std::vector<bool> used(nucleons.size(),false);
  for (size_t i=0;i<nucleons.size();i++) {
    for (size_t j=i+1;j<nucleons.size() && !used[i];j++) {
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
      m_evcounts[std::string("formed_")+(anti ? "dbar_" : "d_")+
                 Deuteron_Coalescence::Name(channel)] += 1.;
      msg_Out()<<"   * "<<(anti ? "anti" : "")<<"deuteron formed via "
               <<Deuteron_Coalescence::Name(channel)<<", p = "<<deut->Momentum()<<"\n";
    }
  }
  for (const auto & c : m_evcounts) {
    m_counts[c.first]  += c.second;
    m_counts2[c.first] += c.second*c.second;
  }
  if (inel0) {
    size_t row = std::min(nV0,mult_nmax-1);
    m_mult_events[row] += 1.;
    m_mult_nch[row]    += double(nch_mid);
    for (const auto & ey : m_ev_yield) {
      for (size_t i=0;i<mult_npt;i++) {
        m_mult_yield[ey.first][row*mult_npt+i]  += ey.second[i];
        m_mult_yield2[ey.first][row*mult_npt+i] += ey.second[i]*ey.second[i];
      }
      m_mult_dndy[ey.first][row] += m_ev_dndy[ey.first];
    }
  }
  return formed ? Return_Value::Success : Return_Value::Nothing;
}

void Deuteron_Formation::CleanUp(const size_t & mode) { m_done = false; }

void Deuteron_Formation::Finish(const std::string &) {
  if (!m_on || m_nevents<=0.) return;
  MakeDir(m_outdir,true);
  // per event and per unit of the histogrammed variable
  for (auto & h : m_histos) {
    h.second->Scale(1./(m_nevents*h.second->BinSize()));
    h.second->Output(m_outdir+"/"+h.first+".dat");
  }
  // pair maps: |p1+p2| (lab) and k (c.m.) bin centres, pairs per event
  double width = map_max/double(map_nbins);
  for (const auto & pm : m_pairmaps) {
    My_Out_File out(m_outdir+"/pairmap_"+pm.first+".dat");
    out.Open();
    *out<<"# |p1+p2| [GeV]  k [GeV]  pairs per event per bin ("
        <<width<<" GeV x "<<width<<" GeV)\n";
    for (size_t ik=0;ik<map_nbins;ik++)
      for (size_t ip=0;ip<map_nbins;ip++)
        *out<<(ip+0.5)*width<<" "<<(ik+0.5)*width<<" "
            <<pm.second[ik*map_nbins+ip]/m_nevents<<"\n";
    out.Close();
  }
  // multiplicity classes: rows = V0M charged multiplicity (last row: overflow),
  // only rows with events; yields not divided by the number of events
  {
    My_Out_File out(m_outdir+"/multiplicity_events.dat");
    out.Open();
    *out<<"# INEL>0 events binned in the V0M charged multiplicity nV0 (V0-A 2.8<=eta<=5.1, "
        <<"V0-C -3.7<=eta<=-1.7; arXiv:1902.09290)\n"
        <<"# nV0  events  sum_dNch(|eta|<0.5)  sum_dN/dy(|y|<0.5): d  dbar  p  pbar\n";
    for (size_t r=0;r<mult_nmax;r++) {
      if (m_mult_events[r]<=0.) continue;
      *out<<r<<" "<<m_mult_events[r]<<" "<<m_mult_nch[r];
      for (std::string sp : {"d","dbar","p","pbar"}) *out<<" "<<m_mult_dndy[sp][r];
      *out<<"\n";
    }
    out.Close();
    double bw = mult_ptmax/double(mult_npt);
    for (std::string sp : {"d","dbar","p","pbar"}) {
      for (std::string sq : {"","_sq"}) {
        const std::vector<double> & vals = sq.empty() ? m_mult_yield[sp] : m_mult_yield2[sp];
        My_Out_File o(m_outdir+"/multiplicity_"+sp+"_invyield_midy"+sq+".dat");
        o.Open();
        *o<<"# rows: nV0 (as multiplicity_events.dat); columns: "<<mult_npt<<" pT bins of "<<bw
          <<" GeV from 0; "<<(sq.empty() ? "sum over events of 1/(2 pi pT) dN/dpT dy, |y|<0.5"
                                          : "sum over events of the squared per-event values")
          <<" [GeV^-2"<<(sq.empty() ? "" : ", squared")<<"]\n";
        double norm = sq.empty() ? bw : bw*bw;
        for (size_t r=0;r<mult_nmax;r++) {
          if (m_mult_events[r]<=0.) continue;
          *o<<r;
          for (size_t i=0;i<mult_npt;i++) *o<<" "<<vals[r*mult_npt+i]/norm;
          *o<<"\n";
        }
        o.Close();
      }
    }
  }
  // per-event yields: key, mean, error (from the spread of per-event yields)
  My_Out_File summary(m_outdir+"/summary.dat");
  summary.Open();
  *summary<<"# Deuteron_Formation: "<<m_name<<", 1/sigma_0 = "
          <<p_deuterons->InvSigma0()<<" 1/b, NSAMP = "<<m_nsamp
          <<", events = "<<m_nevents<<"\n"
          <<"# ALEPH_window: 0.62 < p < 1.03 GeV, |cos(theta)| < 0.95; "
          <<"ALEPH: dbar = (5.9 +- 1.8 +- 0.5)e-6 per hadronic Z event\n"
          <<"# key  per-event-mean  error\n";
  for (const auto & c : m_counts) {
    double mean = c.second/m_nevents;
    double var  = std::max(0.,m_counts2[c.first]/m_nevents-mean*mean);
    *summary<<c.first<<"  "<<mean<<"  "<<sqrt(var/m_nevents)<<"\n";
  }
  summary.Close();
  msg_Info()<<METHOD<<": "<<m_nevents<<" events, results in "<<m_outdir<<"/.\n";
}
