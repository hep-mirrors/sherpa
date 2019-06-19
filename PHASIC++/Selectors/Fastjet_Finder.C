#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"
#include "PHASIC++/Selectors/Selector.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/JadePlugin.hh"

namespace PHASIC {
  class Fastjet_Finder : public Selector_Base {
    double m_ptmin,m_etmin,m_delta_r,m_f,m_eta,m_y;
    int m_nj, m_nb, m_nb2, m_eekt;
    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;
    fastjet::EECambridgePlugin * p_eecamplug;
    fastjet::JadePlugin * p_jadeplug;

  public:
    Fastjet_Finder(Process_Base *const proc,std::string algo,
                   double ptmin, double etmin, double dr, double f,
                   double eta, double y, int nj, int nb, int nb2);

    ~Fastjet_Finder();


    bool   Trigger(ATOOLS::Selector_List &);

    void   BuildCuts(Cut_Data *) {}
  };
}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Fastjet_Finder::Fastjet_Finder(Process_Base *const proc,
                               std::string algo, double ptmin, double etmin,
                               double dr, double f, double eta, double y,
                               int nj, int nb, int nb2) :
  Selector_Base("FastjetFinder",proc), m_ptmin(ptmin), m_etmin(etmin),
  m_delta_r(dr), m_f(f), m_eta(eta), m_y(y),
  m_nj(nj), m_nb(nb), m_nb2(nb2), m_eekt(0), p_jdef(0),
  p_siscplug(NULL), p_eecamplug(NULL), p_jadeplug(NULL)
{
  bool ee(rpa->gen.Beam1().IsLepton() && rpa->gen.Beam2().IsLepton());

  fastjet::RecombinationScheme recom(fastjet::E_scheme);

  size_t sep(algo.find('|'));
  if (sep!=std::string::npos) {
    std::string reco(algo.substr(sep+1));
    algo=algo.substr(0,sep);
    if (reco=="pt") recom=fastjet::pt_scheme;
    if (reco=="pt2") recom=fastjet::pt2_scheme;
    if (reco=="Et") recom=fastjet::Et_scheme;
    if (reco=="Et2") recom=fastjet::Et2_scheme;
    if (reco=="BIpt") recom=fastjet::BIpt_scheme;
    if (reco=="BIpt2") recom=fastjet::BIpt2_scheme;
  }

  fastjet::JetAlgorithm ja(fastjet::kt_algorithm);

  if (algo=="cambridge") ja=fastjet::cambridge_algorithm;
  if (algo=="antikt")    ja=fastjet::antikt_algorithm;
  if (algo=="siscone") p_siscplug=new fastjet::SISConePlugin(m_delta_r,m_f);
  if (ee) {
    if (algo=="eecambridge") p_eecamplug=new fastjet::EECambridgePlugin(dr);
    if (algo=="jade") p_jadeplug=new fastjet::JadePlugin();
  }

  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else if (p_eecamplug) p_jdef=new fastjet::JetDefinition(p_eecamplug);
  else if (p_jadeplug) p_jdef=new fastjet::JetDefinition(p_jadeplug);
  else if (ee) {
    p_jdef=new fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    m_eekt=1;
  }
  else p_jdef=new fastjet::JetDefinition(ja,m_delta_r,recom);

  m_smin       = Max(sqr(m_ptmin),sqr(m_etmin));
}


Fastjet_Finder::~Fastjet_Finder() {
  delete p_jdef;
  if (p_siscplug) delete p_siscplug;
  if (p_eecamplug) delete p_eecamplug;
  if (p_jadeplug) delete p_jadeplug;
}


bool Fastjet_Finder::Trigger(Selector_List &sl)
{
  if (m_nj<1) return true;

  DEBUG_FUNC((p_proc?p_proc->Flavours():Flavour_Vector()));

  std::vector<ATOOLS::Vec4D> p(sl.size());
  for (size_t i(0);i<p.size();++i) p[i]=sl[i].Momentum();
  if (sl[0].Flavour().IsLepton()&&rpa->gen.Beam2().IsHadron()) {
    msg_Debugging()<<METHOD<<"(): Boost to Breit frame {\n";
    Vec4D pp(rpa->gen.PBeam(1)), qq(p[0]-p[2]);
    Poincare cms(pp+qq);
    double Q2(-qq.Abs2()), x(Q2/(2.0*pp*qq)), E(sqrt(Q2)/(2.0*x));
    Vec4D P(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
    Vec4D q(0.0,0.0,0.0,2.0*x*E);
    cms.Boost(pp);
    cms.Boost(qq);
    Poincare zrot(pp,-Vec4D::ZVEC);
    zrot.Rotate(pp);
    zrot.Rotate(qq);
    Poincare breit(P+q);
    breit.BoostBack(pp);
    breit.BoostBack(qq);
    if (!IsEqual(pp,P,1.0e-3) || !IsEqual(qq,q,1.0e-3))
      msg_Error()<<METHOD<<"(): Boost error."<<std::endl;
    for (int i(0);i<p.size();++i) {
      msg_Debugging()<<"  "<<i<<": "<<p[i];
      cms.Boost(p[i]);
      zrot.Rotate(p[i]);
      breit.BoostBack(p[i]);
      msg_Debugging()<<" -> "<<p[i]<<"\n";
    }
    msg_Debugging()<<"}\n";
  }

  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<sl.size();++i) {
    if (ToBeClustered(sl[i].Flavour(), (m_nb>0 || m_nb2>0))) {
      input.push_back(MakePseudoJet(sl[i].Flavour(), p[i]));
    }
  }
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=fastjet::sorted_by_pt(cs.inclusive_jets());
  msg_Debugging()<<"njets(ini)="<<jets.size()<<std::endl;

  if (m_eekt) {
    int n(0);
    for (size_t i(0);i<input.size();++i)
      if (cs.exclusive_dmerge_max(i)>sqr(m_ptmin)) ++n;
    return (1-m_sel_log->Hit(1-(n>=m_nj)));
  }

  int n(0), nb(0), nb2(0);
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    msg_Debugging()<<"Jet "<<i<<": pT="<<pj.PPerp()<<", |eta|="<<dabs(pj.Eta())
                   <<", |y|="<<dabs(pj.Y())<<std::endl;
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin &&
	(m_eta==100 || dabs(pj.Eta())<m_eta) &&
	(m_y==100 || dabs(pj.Y())<m_y)) {
      n++;
      if (BTag(jets[i], 1)) nb++;
      if (BTag(jets[i], 2)) nb2++;
    }
  }
  msg_Debugging()<<"njets(fin)="<<n<<std::endl;

  bool trigger(true);
  if (n<m_nj)   trigger=false;
  if (nb<m_nb)   trigger=false;
  if (nb2<m_nb2) trigger=false;

  if (!trigger) {
    msg_Debugging()<<"Point discarded by jet finder"<<std::endl;
  } else {
    msg_Debugging()<<"Point passed"<<std::endl;
  }

  return (1-m_sel_log->Hit(1-trigger));
}


DECLARE_GETTER(Fastjet_Finder,"FastjetFinder",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Finder>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings["FastjetFinder"];
  const auto algo = s["Algorithm"].SetDefault("").Get<std::string>();
  const auto n = s["N"].SetDefault(0).Get<size_t>();
  const auto ptmin = s["PTMin"].SetDefault(0.0).Get<double>();
  const auto etmin = s["ETMin"].SetDefault(0.0).Get<double>();
  const auto dr = s["DR"].SetDefault(0.4).Get<double>();
  const auto f = s["f"].SetDefault(0.75).Get<double>();
  const auto eta = s["EtaMax"].SetDefault(100.0).Get<double>();
  const auto y = s["YMax"].SetDefault(100.0).Get<double>();
  const auto nb = s["Nb"].SetDefault(-1).Get<int>();
  const auto nb2 = s["Nb2"].SetDefault(-1).Get<int>();
#ifndef USING__FASTJET__3
  if (nb>0 || nb2>0) THROW(fatal_error, "b-tagging needs FastJet >= 3.0.");
#endif
  Fastjet_Finder *jf(new Fastjet_Finder(key.p_proc,
                                        algo,
                                        ptmin, etmin,
                                        dr, f,
                                        eta, y,
                                        n, nb, nb2));
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Finder>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FastjetFinder:\n"
     <<width<<"  Algorithm: kt (default)|antikt|cambridge|siscone (hadron colliders)\n"
     <<width<<"  Algorithm: eekt (default)|jade|eecambridge|siscone (lepton colliders)\n"
     <<width<<"  N: number of jets\n"
     <<width<<"  # optional settings:\n"
     <<width<<"  PTMin: minimum jet pT\n"
     <<width<<"  ETMin: minimum jet eta\n"
     <<width<<"  DR: jet distance parameter\n"
     <<width<<"  f: Siscone f parameter (default: 0.75)\n"
     <<width<<"  EtaMax: maximum jet eta (default: 100)\n"
     <<width<<"  YMax: maximum jet rapidity (default: 100)\n"
     <<width<<"  Nb: number of jets with b quarks\n"
     <<width<<"  Nb2: number of jets with non-vanishing b content";
}

#endif
