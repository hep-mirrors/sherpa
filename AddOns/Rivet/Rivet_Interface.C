#include "SHERPA/Tools/Analysis_Interface.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "SHERPA/Single_Events/Event_Handler.H"

// TODO:
// to be done properly
#define RIVET_ENABLE_HEPMC_3

#ifdef USING__RIVET3
#include "Rivet/Config/RivetConfig.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"
#include "YODA/Config/BuildConfig.h"
#include "YODA/AnalysisObject.h"
#if defined(RIVET_ENABLE_HEPMC_3) || (RIVET_VERSION_CODE >= 30200)
#include "SHERPA/Tools/HepMC3_Interface.H"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenCrossSection.h"
#define SHERPA__HepMC_Interface SHERPA::HepMC3_Interface
#define HEPMCNS HepMC3
#define HEPMC_HAS_CROSS_SECTION
#else
#include "SHERPA/Tools/HepMC2_Interface.H"
#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/WeightContainer.h"
#include "HepMC/HepMCDefs.h"
#define SHERPA__HepMC_Interface SHERPA::HepMC2_Interface
#define HEPMCNS HepMC
#endif

#define USING__Rivet_MPI_Merge

namespace SHERPARIVET {
  typedef std::pair<std::string, int> RivetMapKey;
  typedef std::map<RivetMapKey, Rivet::AnalysisHandler*> Rivet_Map;

  class Rivet_Interface: public SHERPA::Analysis_Interface {
  private:

    std::string m_outpath, m_tag;

    size_t m_nevt;
    bool   m_finished;
    bool   m_splitjetconts, m_splitSH, m_splitpm, m_splitcoreprocs, m_usehepmcshort;

    Rivet_Map         m_rivet;
    SHERPA__HepMC_Interface       m_hepmc2;
    std::vector<ATOOLS::btp::code> m_ignoreblobs;
    std::map<std::string,size_t>   m_weightidxmap;
#if defined(USING__RIVET3) && defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
#ifndef  RIVET_ENABLE_HEPMC_3
    HepMC::GenEvent m_lastevent;
#else
    HepMC3::GenEvent m_lastevent;
#endif
#endif

    Rivet::AnalysisHandler* GetRivet(std::string proc,
                                     int jetcont);
    std::string             GetCoreProc(const std::string& proc);
    std::string             OutputPath(const Rivet_Map::key_type& key);

  public:
    Rivet_Interface(const std::string &outpath,
                    const std::vector<ATOOLS::btp::code> &ignoreblobs,
                    const std::string &tag);
    ~Rivet_Interface();

    bool Init();
    bool Run(ATOOLS::Blob_List *const bl);
    bool Finish();

    void ShowSyntax(const int i);
  };

  class RivetShower_Interface: public Rivet_Interface {};
  class RivetME_Interface: public Rivet_Interface {};
}


using namespace SHERPARIVET;
using namespace SHERPA;
using namespace ATOOLS;
using namespace Rivet;


Rivet_Interface::Rivet_Interface(const std::string &outpath,
                                 const std::vector<btp::code> &ignoreblobs,
                                 const std::string& tag) :
  Analysis_Interface("Rivet"),
  m_outpath(outpath), m_tag(tag),
  m_nevt(0), m_finished(false),
  m_splitjetconts(false), m_splitSH(false), m_splitpm(false),
  m_splitcoreprocs(false),
  m_ignoreblobs(ignoreblobs)
{
  if (m_outpath[m_outpath.size()-1]=='/')
    m_outpath=m_outpath.substr(0,m_outpath.size()-1);
  if (m_outpath.rfind('/')!=std::string::npos)
    MakeDir(m_outpath.substr(0,m_outpath.rfind('/')), true);
#ifdef USING__MPI
  if (mpi->Rank()==0) {
#endif
    if (m_outpath.rfind('/')!=std::string::npos)
      MakeDir(m_outpath.substr(0,m_outpath.rfind('/')));
#ifdef USING__MPI
  }
#if ! ( defined(USING__RIVET3) && defined(USING__MPI) && defined(USING__Rivet_MPI_Merge) )
  if (mpi->Size()>1)
    m_outpath.insert(m_outpath.length(),"_"+rpa->gen.Variable("RNG_SEED"));
#endif
#endif
}

Rivet_Interface::~Rivet_Interface()
{
  if (!m_finished) Finish();
  for (Rivet_Map::iterator it(m_rivet.begin());
       it!=m_rivet.end();++it) {
    delete it->second;
  }
  m_rivet.clear();
}

AnalysisHandler* Rivet_Interface::GetRivet(std::string proc,
                                           int jetcont)
{
  DEBUG_FUNC(proc<<" "<<jetcont);
  RivetMapKey key = std::make_pair(proc, jetcont);
  Rivet_Map::iterator it=m_rivet.find(key);
  if (it==m_rivet.end()) {
    msg_Debugging()<<"create new "<<key.first<<" "<<key.second<<std::endl;
    m_rivet[key] = new AnalysisHandler();
    Scoped_Settings s{ Settings::GetMainSettings()[m_tag] };
    m_rivet[key]->addAnalyses(s["ANALYSES"].SetSynonyms({"ANALYSIS", "-a", "--analyses"})
                       .SetDefault<std::vector<std::string>>({}).GetVector<std::string>());
#if RIVET_VERSION_CODE >= 30200
    m_rivet[key]->setCheckBeams(!s["IGNORE_BEAMS"].SetSynonyms({"IGNOREBEAMS", "--ignore-beams"}).SetDefault(0).Get<int>());
    m_rivet[key]->matchWeightNames(s["MATCH_WEIGHTS"].SetSynonyms({"--match-weights"}).SetDefault("").Get<std::string>());
    m_rivet[key]->unmatchWeightNames(s["UNMATCH_WEIGHTS"].SetSynonyms({"--unmatch-weights"}).SetDefault("^EXTRA__.*,^IRREG__.*").Get<std::string>());
    m_rivet[key]->setFinalizePeriod(OutputPath(key), s["HISTO_INTERVAL"].SetSynonyms({"--histo-interval"}).SetDefault(0).Get<size_t>());
#else
    m_rivet[key]->setIgnoreBeams(s["IGNORE_BEAMS"].SetSynonyms({"IGNOREBEAMS", "--ignore-beams"}).SetDefault(0).Get<int>());
    m_rivet[key]->selectMultiWeights(s["MATCH_WEIGHTS"].SetSynonyms({"--match-weights"}).SetDefault("").Get<std::string>());
    m_rivet[key]->deselectMultiWeights(s["UNMATCH_WEIGHTS"].SetSynonyms({"--unmatch-weights"}).SetDefault("^EXTRA__.*,^IRREG__.*").Get<std::string>());
    m_rivet[key]->setAODump(OutputPath(key), s["HISTO_INTERVAL"].SetSynonyms({"--histo-interval"}).SetDefault(0).Get<size_t>());
#endif
    m_rivet[key]->skipMultiWeights(s["SKIP_WEIGHTS"].SetSynonyms({"SKIPWEIGHTS", "--skip-weights"}).SetDefault(0).Get<int>());
    m_rivet[key]->setNominalWeightName(s["NOMINAL_WEIGHT"].SetSynonyms({"--nominal-weight"}).SetDefault("").Get<std::string>());
    m_rivet[key]->setWeightCap(s["WEIGHT_CAP"].SetSynonyms({"--weight-cap"}).SetDefault(0.0).Get<double>());
    m_rivet[key]->setNLOSmearing(s["NLO_SMEARING"].SetSynonyms({"--nlo-smearing"}).SetDefault(0.0).Get<double>());
    Log::setLevel("Rivet", s["-l"].SetDefault(20).Get<int>());
  }
  return m_rivet[key];
}

std::string Rivet_Interface::GetCoreProc(const std::string& proc)
{
  DEBUG_FUNC(proc);
  size_t idx=5;
  std::vector<ATOOLS::Flavour> flavs;
  while (idx<proc.size()) {
    std::string fl(1, proc[idx]);
    if (fl=="_") {
      ++idx;
      continue;
    }
    for (++idx; idx<proc.size(); ++idx) {
      if (proc[idx]=='_') break;
      fl+=proc[idx];
    }
    bool bar(false);
    if (fl.length()>1) {
      if (fl.back()=='b') {
        fl.pop_back();
        bar=true;
      }
      else if ((fl[0]=='W' || fl[0]=='H')) {
        if (fl.back()=='-') {
          fl.back()='+';
          bar=true;
        }
      }
      else if (fl.back()=='+') {
        fl.back()='-';
        bar=true;
      }
    }
    Flavour flav(s_kftable.KFFromIDName(fl));
    if (bar) flav=flav.Bar();
    flavs.push_back(flav);
  }

  std::vector<Flavour> nojetflavs;
  for (size_t i=2; i<flavs.size(); ++i) {
    if (!Flavour(kf_jet).Includes(flavs[i])) nojetflavs.push_back(flavs[i]);
  }

  std::vector<Flavour> noewjetflavs;
  for (size_t i=0; i<nojetflavs.size(); ++i) {
    if (!Flavour(kf_ewjet).Includes(nojetflavs[i])) noewjetflavs.push_back(nojetflavs[i]);
  }

  std::vector<Flavour> finalflavs;
  // start with initial state
  for (size_t i=0; i<2; ++i) {
    if (Flavour(kf_jet).Includes(flavs[i]))
      finalflavs.push_back(Flavour(kf_jet));
    else if (Flavour(kf_ewjet).Includes(flavs[i]))
      finalflavs.push_back(Flavour(kf_ewjet));
    else
      finalflavs.push_back(flavs[i]);
  }
  // add all non-jet and non-ewjet particles
  for (size_t i=0; i<noewjetflavs.size(); ++i) {
    finalflavs.push_back(noewjetflavs[i]);
  }
  // add all ewjet particles
  for (size_t i=0; i<nojetflavs.size()-noewjetflavs.size(); ++i) {
    if (finalflavs.size()>3) break;
    finalflavs.push_back(Flavour(kf_ewjet));
  }
  // add all jet particles
  for (size_t i=0; i<flavs.size()-2-nojetflavs.size(); ++i) {
    if (finalflavs.size()>3) break;
    finalflavs.push_back(Flavour(kf_jet));
  }

  std::string ret;
  for (size_t i=0; i<finalflavs.size(); ++i) {
    ret+=finalflavs[i].IDName();
    ret+="__";
  }
  while (!ret.empty() && ret.back()=='_') {
    ret.pop_back();
  }

  DEBUG_VAR(ret);
  return ret;
}

bool Rivet_Interface::Init()
{
  if (m_nevt==0) {
    Scoped_Settings s{ Settings::GetMainSettings()[m_tag] };
    m_splitjetconts = s["JETCONTS"].SetDefault(0).Get<int>();
    m_splitSH = s["SPLITSH"].SetDefault(0).Get<int>();
    m_splitpm = s["SPLITPM"].SetDefault(0).Get<int>();
    m_splitcoreprocs = s["SPLITCOREPROCS"].SetDefault(0).Get<int>();
    m_usehepmcshort = s["USE_HEPMC_SHORT"].SetDefault(0).Get<int>();
    if (m_usehepmcshort && m_tag!="RIVET" && m_tag!="RIVETSHOWER") {
      THROW(fatal_error, "Internal error.");
    }

    // configure HepMC interface
    for (size_t i=0; i<m_ignoreblobs.size(); ++i) {
      m_hepmc2.Ignore(m_ignoreblobs[i]);
    }
    m_hepmc2.SetHepMCNamedWeights(
        s["USE_HEPMC_NAMED_WEIGHTS"].SetDefault(true).Get<bool>());
    m_hepmc2.SetHepMCExtendedWeights(
        s["USE_HEPMC_EXTENDED_WEIGHTS"].SetDefault(false).Get<bool>());
    m_hepmc2.SetHepMCTreeLike(
        s["USE_HEPMC_TREE_LIKE"].SetDefault(false).Get<bool>());
  }
  return true;
}

bool Rivet_Interface::Run(ATOOLS::Blob_List *const bl)
{
  DEBUG_FUNC("");
  Particle_List pl=bl->ExtractParticles(1);
  for (Particle_List::iterator it=pl.begin(); it!=pl.end(); ++it) {
    if ((*it)->Momentum().Nan()) {
      msg_Error()<<METHOD<<" encountered NaN in momentum. Ignoring event:"
                 <<std::endl<<*bl<<std::endl;
      return true;
    }
  }

#ifndef  RIVET_ENABLE_HEPMC_3
  HepMC::GenEvent event;
#else
  HepMC3::GenEvent event;
#endif
  if (m_usehepmcshort)  m_hepmc2.Sherpa2ShortHepMC(bl, event);
  else                  m_hepmc2.Sherpa2HepMC(bl, event);
  std::vector<HEPMCNS::GenEvent*> subevents(m_hepmc2.GenSubEventList());
#ifdef  RIVET_ENABLE_HEPMC_3
  m_hepmc2.AddCrossSection(event, p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
#else
  m_hepmc2.AddCrossSection(event, p_eventhandler->TotalNominalXS());
#endif
#if defined(USING__RIVET3) && defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
#ifndef  RIVET_ENABLE_HEPMC_3
  // to be done properly
  // HepMC::GenCrossSection xs;
  // xs.set_cross_section(p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
  if (m_lastevent.vertices_begin()==m_lastevent.vertices_end()) {
#else
  // to be done properly
  // std::shared_ptr<HepMC3::GenCrossSection> xs=std::make_shared<HepMC3::GenCrossSection>();
  // xs->set_cross_section(p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
  if (m_lastevent.vertices().empty()) {
#endif
    m_lastevent=event;
    for (size_t i(0);i<m_lastevent.weights().size();++i) m_lastevent.weights()[i]=0;
  }
  // to be done properly
  // m_lastevent.set_cross_section(xs);
#endif

  if (subevents.size()) {
    for (size_t i(0);i<subevents.size();++i) {
      GetRivet("",0)->analyze(*subevents[i]);
    }
    m_hepmc2.DeleteGenSubEventList();
  }
  else {
    GetRivet("",0)->analyze(event);
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    size_t parts=0;
    if (sp) {
      std::string multi(sp?sp->TypeSpec():"");
      if (multi[3]=='_') multi=multi.substr(2,1);
      else multi=multi.substr(2,2);
      parts=ToType<size_t>(multi);
    }
    if (m_splitjetconts && sp) {
      GetRivet("",parts)->analyze(event);
    }
    if (m_splitcoreprocs && sp) {
      GetRivet(GetCoreProc(sp->TypeSpec()),0)->analyze(event);
      if (m_splitjetconts) {
        GetRivet(GetCoreProc(sp->TypeSpec()),parts)->analyze(event);
      }
    }
    if (m_splitSH && sp) {
      std::string typespec=sp->TypeSpec();
      typespec=typespec.substr(typespec.length()-2, 2);
      std::string type="";
      if (typespec=="+S") type="S";
      else if (typespec=="+H") type="H";
      
      if (type!="") {
        GetRivet(type,0)->analyze(event);
        if (m_splitjetconts) {
          GetRivet(type,parts)->analyze(event);
        }
      }
    }
    if (m_splitpm) {
      GetRivet(event.weights()[0]<0?"M":"P",0)->analyze(event);
    }
  }

  ++m_nevt;
  return true;
}

std::string Rivet_Interface::OutputPath(const Rivet_Map::key_type& key)
{
  std::string out = m_outpath;
  if (key.first!="") out+="."+key.first;
  if (key.second!=0) out+=".j"+ToString(key.second);
  out+=".yoda";
#ifdef HAVE_LIBZ
  out+=".gz";
#endif
  return out;
}

bool Rivet_Interface::Finish()
{
  PRINT_FUNC(m_outpath);
  #if RIVET_VERSION_CODE >= 30200
  GetRivet("",0)->pushToPersistent();
  #else
  GetRivet("",0)->finalize();
  #endif
  const double nomsumw = GetRivet("",0)->sumW();
  #ifdef  RIVET_ENABLE_HEPMC_3
  const double nomxsec = p_eventhandler->TotalXS().Nominal();
  const double nomxerr = p_eventhandler->TotalErr().Nominal();
  #else
  const double nomxsec = p_eventhandler->TotalNominalXS().value;
  const double nomxerr = 0.0;
  #endif

  // additional Rivet instances are used e.g. for splits
  // into H+S events or jet multiplicities -> these only
  // get to see a subset of the events and need to be
  // re-scaled to cross-section of the complete run
  const bool needs_rescaling = m_rivet.size() > 1;
  for (auto& it : m_rivet) {
    if (needs_rescaling) {
      // first collapse the event group,
      // then scale the cross-section
      // before finalizing
      #if RIVET_VERSION_CODE >= 30200
      it.second->pushToPersistent();
      #else
      it.second->finalize();
      #endif

      // determine the weight fraction seen by this Rivet run
      const double wgtfrac = it.second->sumW()/nomsumw;
      // rescale nominal cross-section
      const double thisxs  = nomxsec*wgtfrac;
      const double thiserr = nomxerr*wgtfrac;
      it.second->setCrossSection(thisxs, thiserr);
    }
    it.second->finalize();
    it.second->writeData(OutputPath(it.first));
  }
  m_finished=true;
  return true;
}

void Rivet_Interface::ShowSyntax(const int i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
    <<"   RIVET: {\n\n"
    <<"     --analyses: [<ana_1>, <ana_2>]  # analyses to run\n"
    <<"     # Optional parameters: Please refer to manual\n"
    <<"}"<<std::endl;
}

DECLARE_GETTER(Rivet_Interface,"Rivet",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,Rivet_Interface>::
operator()(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath.back()=='/') {
    outpath.pop_back();
  }
  std::vector<btp::code> ignoreblobs;
  ignoreblobs.push_back(btp::Unspecified);
  return new Rivet_Interface(outpath,
                             ignoreblobs,
                             "RIVET");
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,Rivet_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface";
}


DECLARE_GETTER(RivetShower_Interface,"RivetShower",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,RivetShower_Interface>::
operator()(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath.back()=='/') {
    outpath.pop_back();
  }
  std::vector<btp::code> ignoreblobs;
  ignoreblobs.push_back(btp::Unspecified);
  ignoreblobs.push_back(btp::Fragmentation);
  ignoreblobs.push_back(btp::Hadron_Decay);
  ignoreblobs.push_back(btp::Hadron_Mixing);
  return new Rivet_Interface(outpath + ".SL",
                             ignoreblobs,
                             "RIVETSHOWER");
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,RivetShower_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface on top of shower level events.";
}


DECLARE_GETTER(RivetME_Interface,"RivetME",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,RivetME_Interface>::
operator()(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath.back()=='/') {
    outpath.pop_back();
  }
  std::vector<btp::code> ignoreblobs;
  ignoreblobs.push_back(btp::Unspecified);
  ignoreblobs.push_back(btp::Fragmentation);
  ignoreblobs.push_back(btp::Hadron_Decay);
  ignoreblobs.push_back(btp::Hadron_Mixing);
  ignoreblobs.push_back(btp::Shower);
  ignoreblobs.push_back(btp::Hadron_To_Parton);
  ignoreblobs.push_back(btp::Hard_Collision);
  ignoreblobs.push_back(btp::QED_Radiation);
  ignoreblobs.push_back(btp::Soft_Collision);
  return new Rivet_Interface(outpath + ".ME",
                             ignoreblobs,
                             "RIVETME");
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,RivetME_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface on top of ME level events.";
}

#endif
// end of Rivet_Interface for Rivet3











#ifdef USING__RIVET2
#include "SHERPA/Tools/HepMC2_Interface.H"
#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/WeightContainer.h"
#include "HepMC/HepMCDefs.h"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"

namespace SHERPARIVET {
  typedef std::pair<std::string, int> RivetMapKey;
  typedef std::map<RivetMapKey, Rivet::AnalysisHandler*> Rivet_Map;

  class Rivet_Scale_Variation {
  private:
    std::string m_name;
    Rivet_Map   m_rivetmap;
    double      m_wgt,m_n,m_sum,m_sum2,m_tempn,m_tempsum;
    std::vector<double> m_rswgts;
  public:
    Rivet_Scale_Variation(std::string name="");
    ~Rivet_Scale_Variation();

    void   AddPoint(const double& wgt, const double& n, size_t xsmode=0);
    void   SynchroniseCrossSection();
    double TotalXS()  const;
    double TotalVar() const;
    double TotalErr() const;

    inline Rivet_Map&  RivetMap() { return m_rivetmap; }

    inline std::string Name()         const { return m_name; }
    inline double      Weight()       const { return m_wgt; }
    inline double      Weight(size_t i) const { return m_rswgts[i]; }
    inline double      SumOfWeights() const { return m_sum; }

    inline void ResetRSWeights() { m_rswgts.clear(); }
  };

  typedef std::map<std::string, Rivet_Scale_Variation *> RivetScaleVariationMap;

  class Rivet_Interface: public SHERPA::Analysis_Interface {
  private:

    std::string m_outpath, m_tag;
    std::vector<std::string> m_analyses;

    size_t m_nevt;
    bool   m_finished;
    bool   m_splitjetconts, m_splitSH, m_splitpm,
           m_splitcoreprocs, m_splitvariations,
           m_ignorebeams, m_usehepmcshort,
           m_printsummary,
           m_evtbyevtxs;
    size_t m_hepmcoutputprecision, m_xsoutputprecision;

    RivetScaleVariationMap         m_rivet;
    SHERPA::HepMC2_Interface       m_hepmc2;
    std::vector<ATOOLS::btp::code> m_ignoreblobs;
    std::map<std::string,size_t>   m_weightidxmap;

    void ExtractVariations(const HepMC::GenEvent& evt,
                           const std::vector<HepMC::GenEvent*>& subevents);
    void ExtractVariations(const HepMC::GenEvent& evt);
    void SetEventWeight(const Rivet_Scale_Variation* rsv,
                        HepMC::GenEvent& evt,
                        const int& idx=-1);
    void ResetRivetScaleVariationMapRSWeights();

    Rivet::AnalysisHandler* GetRivet(Rivet_Map& rm, std::string proc,
                                     int jetcont);
    std::string             GetCoreProc(const std::string& proc);

  public:
    Rivet_Interface(const std::string &outpath,
                    const std::vector<ATOOLS::btp::code> &ignoreblobs,
                    const std::string &tag);
    ~Rivet_Interface();

    bool Init();
    bool Run(ATOOLS::Blob_List *const bl);
    bool Finish();

    void ShowSyntax(const int i);
  };

  class RivetShower_Interface: public Rivet_Interface {};
  class RivetME_Interface: public Rivet_Interface {};
}


using namespace SHERPARIVET;
using namespace SHERPA;
using namespace ATOOLS;
using namespace Rivet;


Rivet_Scale_Variation::Rivet_Scale_Variation(std::string name) :
  m_name(name), m_rivetmap(),
  m_wgt(0.), m_n(0.), m_sum(0.), m_sum2(0.), m_tempn(0.), m_tempsum(0.),
  m_rswgts(0,0.)
{
}

Rivet_Scale_Variation::~Rivet_Scale_Variation()
{
  for (Rivet_Map::iterator it=m_rivetmap.begin();it!=m_rivetmap.end();++it) {
    delete it->second;
  }
}

void Rivet_Scale_Variation::AddPoint(const double& wgt, const double& n,
                                     size_t xsmode)
{
  DEBUG_FUNC("wgt="<<wgt<<", n="<<n<<", mode="<<xsmode);
  if      (xsmode==0) {
    m_wgt=wgt;
    m_n+=n;
    m_sum+=wgt;
    m_sum2+=ATOOLS::sqr(wgt);
  }
  else if (xsmode==1) {
    if (m_tempn>0 && m_tempn!=n) THROW(fatal_error,"Inconsistent ntrial.");
    m_rswgts.push_back(wgt);
    m_tempn=n;
    m_tempsum+=wgt;
  }
  else THROW(fatal_error,"Unknown xs-mode.");
}

void Rivet_Scale_Variation::SynchroniseCrossSection()
{
  m_n+=m_tempn;
  m_sum+=m_tempsum;
  m_sum2+=ATOOLS::sqr(m_tempsum);
  m_tempn=m_tempsum=0.;
}

double Rivet_Scale_Variation::TotalXS() const
{
  if (m_n==0.) return 0.;
  return m_sum/m_n;
}

double Rivet_Scale_Variation::TotalVar() const
{
  if (m_n<=1.) return ATOOLS::sqr(TotalXS());
  return (m_sum2-m_sum*m_sum/m_n)/(m_n-1.);
}

double Rivet_Scale_Variation::TotalErr() const
{
  if (m_n<=1.) return TotalXS();
  if (ATOOLS::IsEqual
      (m_sum2*m_n,m_sum*m_sum,1.0e-6)) return 0.;
  return sqrt((m_sum2-m_sum*m_sum/m_n)/(m_n-1.)/m_n);
}

Rivet_Interface::Rivet_Interface(const std::string &outpath,
                                 const std::vector<btp::code> &ignoreblobs,
                                 const std::string& tag) :
  Analysis_Interface("Rivet"),
  m_outpath(outpath), m_tag(tag),
  m_nevt(0), m_finished(false),
  m_splitjetconts(false), m_splitSH(false), m_splitpm(false),
  m_splitcoreprocs(false), m_splitvariations(true),
  m_ignoreblobs(ignoreblobs),
  m_printsummary(true), m_evtbyevtxs(false),
  m_hepmcoutputprecision(15), m_xsoutputprecision(6)
{
  if (m_outpath[m_outpath.size()-1]=='/')
    m_outpath=m_outpath.substr(0,m_outpath.size()-1);
  if (m_outpath.rfind('/')!=std::string::npos)
    MakeDir(m_outpath.substr(0,m_outpath.rfind('/')));
#ifdef USING__MPI
  if (mpi->Size()>1) {
    m_outpath.insert(m_outpath.length(),"_"+rpa->gen.Variable("RNG_SEED"));
  }
#endif
}

Rivet_Interface::~Rivet_Interface()
{
  if (!m_finished) Finish();
  for (RivetScaleVariationMap::iterator it(m_rivet.begin());
       it!=m_rivet.end();++it) {
    delete it->second;
  }
  m_rivet.clear();
}

void Rivet_Interface::ExtractVariations
(const HepMC::GenEvent& evt,const std::vector<HepMC::GenEvent*>& subevents)
{
  DEBUG_FUNC("# of subevts: "<<subevents.size());
  if (subevents.size()) {
    for (size_t i(0);i<subevents.size();++i) {
      ExtractVariations(*subevents[i]);
    }
  }
  else ExtractVariations(evt);
}

void Rivet_Interface::ExtractVariations(const HepMC::GenEvent& evt)
{
  DEBUG_FUNC("");
  const HepMC::WeightContainer& wc(evt.weights());
  std::map<std::string,double> wgtmap;
  double ntrials(1.);
  size_t xstype(0);
  // lookup all evt-wgts with a variation weight name
  // at the moment the only way to do that is to filter the printout
  // accuracy limited to print out accu of 6 digits, must suffice
  MyStrStream str;
  str.precision(m_hepmcoutputprecision);
  wc.print(str);

  // need a temp object first, as we need to get ntrials first
  while (str) {
    double wgt(0.);
    std::string cur("");
    str>>cur;
    if (cur.length()==0) continue;
    // weight is between "," and trailing bracket
    // name is between leading bracket and ","
    const auto wgtstart = cur.find(",") + 1;
    const auto wgtend = cur.find(")") - 1;
    wgt = ToType<double>(cur.substr(wgtstart, wgtend - wgtstart + 1));
    cur = cur.substr(1, wgtstart - 2);
    if (m_splitvariations && m_hepmc2.StartsLikeVariationName(cur)) {
      wgtmap[cur]=wgt;
    }
    else if (cur=="Weight")  wgtmap["nominal"]=wgt;
    else if (cur=="EXTRA__NTrials") ntrials=wgt;
    else if (cur=="IRREG__Reweight_Type" && ((int)wgt)&64) xstype=1;
  }
  wgtmap["nominal"]=wc[0];
  ntrials=wc[3];
  xstype=(((wc.size()==5&&((int)wc[4]&64))||(wc.size()==11&&((int)wc[10]&64)))?1:0);
  if (msg_LevelIsDebugging()) {
    for (std::map<std::string,double>::iterator wit(wgtmap.begin());
         wit!=wgtmap.end();++wit)
      msg_Out()<<wit->first<<" : "<<wit->second<<std::endl;
  }
  // now construct or fill into the scale variation map
  for (std::map<std::string,double>::iterator wit(wgtmap.begin());
       wit!=wgtmap.end();++wit) {
    RivetScaleVariationMap::iterator rit=m_rivet.find(wit->first);
    if (rit==m_rivet.end()) {
      msg_Debugging()<<"creating new entry in m_rivet"<<std::endl;
      m_rivet[wit->first]=new Rivet_Scale_Variation(wit->first);
      m_rivet[wit->first]->AddPoint(wit->second,ntrials,xstype);
    }
    else rit->second->AddPoint(wit->second,ntrials,xstype);
  }
  if (msg_LevelIsDebugging()) {
    for (RivetScaleVariationMap::iterator rit(m_rivet.begin());
         rit!=m_rivet.end();++rit)
      msg_Out()<<rit->first<<" : "<<rit->second->Weight()<<std::endl;
  }
}

void Rivet_Interface::SetEventWeight(const Rivet_Scale_Variation* rsv,
                                     HepMC::GenEvent& evt, const int& idx)
{
  double wgt(idx<0?rsv->Weight():rsv->Weight(idx));
  DEBUG_FUNC(rsv->Name()<<": "<<wgt);
  evt.weights()[0]=wgt;
  evt.cross_section()->set_cross_section(rsv->TotalXS(),rsv->TotalErr());
}

void Rivet_Interface::ResetRivetScaleVariationMapRSWeights()
{
  for (RivetScaleVariationMap::iterator rit(m_rivet.begin());
       rit!=m_rivet.end();++rit) rit->second->ResetRSWeights();
}

AnalysisHandler* Rivet_Interface::GetRivet(Rivet_Map& rm, std::string proc,
                                           int jetcont)
{
  DEBUG_FUNC(proc<<" "<<jetcont);
  RivetMapKey key = std::make_pair(proc, jetcont);
  Rivet_Map::iterator it=rm.find(key);
  if (it!=rm.end()) {
    msg_Debugging()<<"found "<<key.first<<" "<<key.second<<std::endl;
    return it->second;
  }
  else {
    msg_Debugging()<<"create new "<<key.first<<" "<<key.second<<std::endl;
    AnalysisHandler* rivet(new AnalysisHandler());
    rivet->setIgnoreBeams(m_ignorebeams);
    rivet->addAnalyses(m_analyses);
    rm.insert(std::make_pair(key, rivet));
    msg_Debugging()<<"now "<<rm.size()<<" in "
                   <<key.first<<" "<<key.second<<std::endl;
    return rivet;
  }
}

std::string Rivet_Interface::GetCoreProc(const std::string& proc)
{
  DEBUG_FUNC(proc);
  size_t idx=5;
  std::vector<ATOOLS::Flavour> flavs;
  while (idx<proc.size()) {
    std::string fl(1, proc[idx]);
    if (fl=="_") {
      ++idx;
      continue;
    }
    for (++idx; idx<proc.size(); ++idx) {
      if (proc[idx]=='_') break;
      fl+=proc[idx];
    }
    bool bar(false);
    if (fl.length()>1) {
      if (fl.back()=='b') {
        fl.pop_back();
        bar=true;
      }
      else if ((fl[0]=='W' || fl[0]=='H')) {
        if (fl.back()=='-') {
          fl.back()='+';
          bar=true;
        }
      }
      else if (fl.back()=='+') {
        fl.back()='-';
        bar=true;
      }
    }
    Flavour flav(s_kftable.KFFromIDName(fl));
    if (bar) flav=flav.Bar();
    flavs.push_back(flav);
  }

  std::vector<Flavour> nojetflavs;
  for (size_t i=2; i<flavs.size(); ++i) {
    if (!Flavour(kf_jet).Includes(flavs[i])) nojetflavs.push_back(flavs[i]);
  }

  std::vector<Flavour> noewjetflavs;
  for (size_t i=0; i<nojetflavs.size(); ++i) {
    if (!Flavour(kf_ewjet).Includes(nojetflavs[i])) noewjetflavs.push_back(nojetflavs[i]);
  }

  std::vector<Flavour> finalflavs;
  // start with initial state
  for (size_t i=0; i<2; ++i) {
    if (Flavour(kf_jet).Includes(flavs[i]))
      finalflavs.push_back(Flavour(kf_jet));
    else if (Flavour(kf_ewjet).Includes(flavs[i]))
      finalflavs.push_back(Flavour(kf_ewjet));
    else
      finalflavs.push_back(flavs[i]);
  }
  // add all non-jet and non-ewjet particles
  for (size_t i=0; i<noewjetflavs.size(); ++i) {
    finalflavs.push_back(noewjetflavs[i]);
  }
  // add all ewjet particles
  for (size_t i=0; i<nojetflavs.size()-noewjetflavs.size(); ++i) {
    if (finalflavs.size()>3) break;
    finalflavs.push_back(Flavour(kf_ewjet));
  }
  // add all jet particles
  for (size_t i=0; i<flavs.size()-2-nojetflavs.size(); ++i) {
    if (finalflavs.size()>3) break;
    finalflavs.push_back(Flavour(kf_jet));
  }

  std::string ret;
  for (size_t i=0; i<finalflavs.size(); ++i) {
    ret+=finalflavs[i].IDName();
    ret+="__";
  }
  while (!ret.empty() && ret.back()=='_') {
    ret.pop_back();
  }

  DEBUG_VAR(ret);
  return ret;
}

bool Rivet_Interface::Init()
{
  if (m_nevt==0) {
    Scoped_Settings s{ Settings::GetMainSettings()[m_tag] };
    m_splitjetconts = s["JETCONTS"].SetDefault(0).Get<int>();
    m_splitSH = s["SPLITSH"].SetDefault(0).Get<int>();
    m_splitcoreprocs = s["SPLITCOREPROCS"].SetDefault(0).Get<int>();
    m_splitvariations = s["SPLITVARIATIONS"].SetDefault(1).Get<int>();
    m_usehepmcshort = s["USE_HEPMC_SHORT"].SetDefault(0).Get<int>();
    if (m_usehepmcshort && m_tag!="RIVET" && m_tag!="RIVETSHOWER") {
      THROW(fatal_error, "Internal error.");
    }
    m_printsummary = s["PRINT_SUMMARY"].SetDefault(1).Get<int>();
    m_evtbyevtxs = s["EVENTBYEVENTXS"].SetDefault(0).Get<int>();
    m_ignorebeams = s["IGNORE_BEAMS"].SetSynonyms({"IGNOREBEAMS", "--ignore-beams"}).SetDefault(0).Get<int>();

    m_hepmcoutputprecision = s["HEPMC_OUTPUT_PRECISION"].SetDefault(15).Get<int>();
    m_xsoutputprecision = s["XS_OUTPUT_PRECISION"].SetDefault(6).Get<int>();
    Log::setLevel("Rivet", s["-l"].SetDefault(20).Get<int>());
    m_analyses = s["ANALYSES"]
      .SetDefault<std::vector<std::string>>({})
      .SetSynonyms({"ANALYSIS", "-a", "--analyses"})
      .GetVector<std::string>();
    for (size_t i(0);i<m_analyses.size();++i) {
      if (m_analyses[i]==std::string("MC_XS")) break;
      if (i==m_analyses.size()-1) m_analyses.push_back(std::string("MC_XS"));
    }

    // configure HepMC interface
    for (size_t i=0; i<m_ignoreblobs.size(); ++i) {
      m_hepmc2.Ignore(m_ignoreblobs[i]);
    }
    m_hepmc2.SetHepMCNamedWeights(
        s["USE_HEPMC_NAMED_WEIGHTS"].SetDefault(true).Get<bool>());
    m_hepmc2.SetHepMCExtendedWeights(
        s["USE_HEPMC_EXTENDED_WEIGHTS"].SetDefault(false).Get<bool>());
    m_hepmc2.SetHepMCTreeLike(
        s["USE_HEPMC_TREE_LIKE"].SetDefault(false).Get<bool>());
  }
  return true;
}

bool Rivet_Interface::Run(ATOOLS::Blob_List *const bl)
{
  DEBUG_FUNC("");
  Particle_List pl=bl->ExtractParticles(1);
  for (Particle_List::iterator it=pl.begin(); it!=pl.end(); ++it) {
    if ((*it)->Momentum().Nan()) {
      msg_Error()<<METHOD<<" encountered NaN in momentum. Ignoring event:"
                 <<endl<<*bl<<endl;
      return true;
    }
  }

  HepMC::GenEvent event;
  if (m_usehepmcshort)  m_hepmc2.Sherpa2ShortHepMC(bl, event);
  else                  m_hepmc2.Sherpa2HepMC(bl, event);
  std::vector<HepMC::GenEvent*> subevents(m_hepmc2.GenSubEventList());
  // leave this, although will be overwritten later
  m_hepmc2.AddCrossSection(event, p_eventhandler->TotalNominalXS());

  // 1st event build index map, thereafter only lookup
  ExtractVariations(event,subevents);
  for (RivetScaleVariationMap::iterator it(m_rivet.begin());
       it!=m_rivet.end();++it) {
    msg_Debugging()<<"Running rivet for "<<it->first<<" with "
                   <<it->second->RivetMap().size()<<" histograms."<<std::endl;
    Rivet_Map& rivetmap(it->second->RivetMap());
    if (subevents.size()) {
      it->second->SynchroniseCrossSection();
      for (size_t i(0);i<subevents.size();++i) {
        SetEventWeight(it->second,*subevents[i],i);
        GetRivet(rivetmap,"", 0)->analyze(*subevents[i]);
      }
    }
    else {
      SetEventWeight(it->second,event);
      GetRivet(rivetmap,"",0)->analyze(event);
      Blob *sp(bl->FindFirst(btp::Signal_Process));
      size_t parts=0;
      if (sp) {
        std::string multi(sp?sp->TypeSpec():"");
        if (multi[3]=='_') multi=multi.substr(2,1);
        else multi=multi.substr(2,2);
        parts=ToType<size_t>(multi);
      }
      if (m_splitjetconts && sp) {
        GetRivet(rivetmap,"",parts)->analyze(event);
      }
      if (m_splitcoreprocs && sp) {
        GetRivet(rivetmap,GetCoreProc(sp->TypeSpec()),0)
            ->analyze(event);
        if (m_splitjetconts) {
          GetRivet(rivetmap,GetCoreProc(sp->TypeSpec()),parts)
              ->analyze(event);
        }
      }
      if (m_splitSH && sp) {
        std::string typespec=sp->TypeSpec();
        typespec=typespec.substr(typespec.length()-2, 2);
        std::string type="";
        if (typespec=="+S") type="S";
        else if (typespec=="+H") type="H";

        if (type!="") {
          GetRivet(rivetmap,type,0)->analyze(event);
          if (m_splitjetconts) {
            GetRivet(rivetmap,type,parts)->analyze(event);
          }
        }
      }
      if (m_splitpm) {
        GetRivet(rivetmap,event.weights()[0]<0?"M":"P",0)->analyze(event);
      }
    }
  }
  if (subevents.size()) {
    ResetRivetScaleVariationMapRSWeights();
    m_hepmc2.DeleteGenSubEventList();
  }

  ++m_nevt;
  for (RivetScaleVariationMap::iterator it(m_rivet.begin());
       it!=m_rivet.end();++it) {
    msg_Debugging()<<"Checking rivet for "<<it->first<<" with "
                   <<it->second->RivetMap().size()<<" histograms."<<std::endl;
  }
  if (m_evtbyevtxs) {
    for (RivetScaleVariationMap::iterator mit(m_rivet.begin());
         mit!=m_rivet.end();++mit) {
      std::string out=m_outpath;
      if (mit->first!="" && mit->first!="nominal") out += "."+mit->first;
      std::string namestr(m_rivet.size()>1?" for "+mit->first:"");
      std::string output(std::string("**  Total XS")+namestr
                         +std::string(" = ( ")
                         +ToString(mit->second->TotalXS(),m_xsoutputprecision)
                         +std::string(" +- ")
                         +ToString(mit->second->TotalErr(),m_xsoutputprecision)
                         +std::string(" ) pb **"));
      std::string astline(output.size(),'*');
      msg_Info()<<astline<<"\n"<<output<<"\n"<<astline<<std::endl;
    }
  }
  return true;
}

bool Rivet_Interface::Finish()
{
  std::string ending("");
  ending=".yoda";
  for (RivetScaleVariationMap::iterator mit(m_rivet.begin());
       mit!=m_rivet.end();++mit) {
    std::string out=m_outpath;
    if (mit->first!="" && mit->first!="nominal") out += "."+mit->first;
    PRINT_FUNC(out+ending);
    if (m_printsummary) {
      std::string namestr(m_rivet.size()>1?" for "+mit->first:"");
      std::string output(std::string("**  Total XS")+namestr
                         +std::string(" = ( ")
                         +ToString(mit->second->TotalXS(),m_xsoutputprecision)
                         +std::string(" +- ")
                         +ToString(mit->second->TotalErr(),m_xsoutputprecision)
                         +std::string(" ) pb **"));
      std::string astline(output.size(),'*');
      msg_Info()<<astline<<"\n"<<output<<"\n"<<astline<<std::endl;
    }
#if defined(USING__RIVET3) && defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
    // to be done properly
    // const double nomxsec = p_eventhandler->TotalXSMPI();
    // const double nomxerr = p_eventhandler->TotalErrMPI();
#else
    // to be done properly
    // const double nomxsec = p_eventhandler->TotalXS();
    // const double nomxerr = p_eventhandler->TotalErr();
#endif
#if defined(USING__RIVET3) && defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
    // synchronize analyses among MPI processes
    std::string mynames;
    for (Rivet_Map::iterator it=mit->second->RivetMap().begin();
	 it!=mit->second.RivetMap().end();++it) {
      std::string out;
      if (it->first.first!="") out+="."+it->first.first;
      if (it->first.second!=0) out+=".j"+ToString(it->first.second);
      mynames+=out+"|";
    }
    int len(mynames.length()+1);
    mpi->Allreduce(&len,1,MPI_INT,MPI_MAX);
    std::string allnames;
    mynames.reserve(len);
    allnames.reserve(len*mpi->Size()+1);
    mpi->Allgather(&mynames[0],len,MPI_CHAR,&allnames[0],len,MPI_CHAR);
    char *catname = new char[len+1];
    for (size_t i(0);i<mpi->Size();++i) {
      sprintf(catname,"%s",&allnames[len*i]);
      std::string curname(catname);
      for (size_t epos(curname.find('|'));
	   epos<curname.length();epos=curname.find('|')) {
	std::string cur(curname.substr(0,epos)), proc, jets;
	curname=curname.substr(epos+1,curname.length()-epos-1);
	size_t dpos(cur.find('.'));
	if (dpos<cur.length()) {
	  proc=cur.substr(dpos+1,cur.length()-dpos-1);
	  cur=cur.substr(0,dpos);
	  size_t jpos(proc.find(".j"));
	  if (jpos<proc.length()) {
	    jets=proc.substr(jpos+2,proc.length()-jpos-1);
	    proc=proc.substr(0,jpos);
	  }
	  else if (proc[0]=='j' && proc.length()>1) {
	    bool isnumber(true);
	    for (size_t j(1);j<proc.length();++j)
	      if (!isdigit(proc[j])) isnumber=false;
	    if (isnumber) {
	      jets=proc.substr(1,proc.length()-1);
	      proc="";
	    }
	  }
	}
	if (jets=="") jets="0";
	RivetMapKey key = std::make_pair(proc,ToType<int>(jets));
	Rivet_Map::iterator it=m_rivet.find(key);
	if (it==m_rivet.end()) {
	  AnalysisHandler* rivet(new AnalysisHandler());
	  rivet->setIgnoreBeams(m_ignorebeams);
	  rivet->skipMultiWeights(m_skipweights);
	  rivet->addAnalyses(m_analyses);
	  rivet->init(m_lastevent);
	  m_rivet.insert(std::make_pair(key, rivet));
	}
      }
    }
    delete [] catname;
#endif
    for (Rivet_Map::iterator it=mit->second->RivetMap().begin();
         it!=mit->second->RivetMap().end(); ++it) {
      const double wgtfrac = it->second->sumOfWeights()/mit->second->SumOfWeights();
      const double totalxs = it->second->crossSection();
      const double thisxs  = totalxs*wgtfrac;
      it->second->setCrossSection(thisxs);
      std::string jout=out;
      if (it->first.first!="") jout+="."+it->first.first;
      if (it->first.second!=0) jout+=".j"+ToString(it->first.second);
      it->second->finalize();
#if defined(USING__RIVET3) && defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
      std::vector<double> data;
      it->second->serialize(data);
      mpi->Allreduce(&data[0],data.size(),MPI_DOUBLE,MPI_SUM);
      size_t i(0);
      it->second->deserialize(data,i);
      if (i!=data.size()) THROW(fatal_error,"serialization error");
      if (mpi->Rank()==0)
#endif
      it->second->writeData(jout+ending);
    }
  }
  m_finished=true;
  return true;
}

void Rivet_Interface::ShowSyntax(const int i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
    <<"   RIVET: {\n\n"
    <<"     ANALYSES: [<ana_1>, <ana_2>]  # analyses to run\n"
    <<"     # optional parameters:\n"
    <<"     JETCONTS: <0|1>      # perform additional separate analyses for \n"
    <<"                          # each matrix element multiplicity\n"
    <<"     SPLITCOREPROCS: <0|1> # perform additional separate analyses for \n"
    <<"                          # each different core process\n"
    <<"     SPLITSH: <0|1>       # perform additional separate analyses for \n"
    <<"                          # S-MC@NLO S- and H- events\n"
    <<"     IGNOREBEAMS: <0|1>   # tell Rivet to ignore beam information\n"
    <<"     USE_HEPMC_SHORT: <0|1> # use shortened HepMC event format\n"
    <<"     USE_HEPMC_NAMED_WEIGHTS: <true|false> # use named HepMC weights,\n"
    <<"                          # mandatory for scale variations\n"
    <<"     PRINT_SUMMARY: <1|0> # print cross section summary at the end\n"
    <<"     EVENTBYEVENTXS: <0|1> # print cross section event-by-event\n"
    <<"}"<<std::endl;
}

DECLARE_GETTER(Rivet_Interface,"Rivet",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,Rivet_Interface>::
operator()(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath.back()=='/') {
    outpath.pop_back();
  }
  std::vector<btp::code> ignoreblobs;
  ignoreblobs.push_back(btp::Unspecified);
  return new Rivet_Interface(outpath,
                             ignoreblobs,
                             "RIVET");
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,Rivet_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface";
}


DECLARE_GETTER(RivetShower_Interface,"RivetShower",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,RivetShower_Interface>::
operator()(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath.back()=='/') {
    outpath.pop_back();
  }
  std::vector<btp::code> ignoreblobs;
  ignoreblobs.push_back(btp::Unspecified);
  ignoreblobs.push_back(btp::Fragmentation);
  ignoreblobs.push_back(btp::Hadron_Decay);
  ignoreblobs.push_back(btp::Hadron_Mixing);
  return new Rivet_Interface(outpath + ".SL",
                             ignoreblobs,
                             "RIVETSHOWER");
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,RivetShower_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface on top of shower level events.";
}


DECLARE_GETTER(RivetME_Interface,"RivetME",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,RivetME_Interface>::
operator()(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath.back()=='/') {
    outpath.pop_back();
  }
  std::vector<btp::code> ignoreblobs;
  ignoreblobs.push_back(btp::Unspecified);
  ignoreblobs.push_back(btp::Fragmentation);
  ignoreblobs.push_back(btp::Hadron_Decay);
  ignoreblobs.push_back(btp::Hadron_Mixing);
  ignoreblobs.push_back(btp::Shower);
  ignoreblobs.push_back(btp::Hadron_To_Parton);
  ignoreblobs.push_back(btp::Hard_Collision);
  ignoreblobs.push_back(btp::QED_Radiation);
  ignoreblobs.push_back(btp::Soft_Collision);
  return new Rivet_Interface(outpath + ".ME",
                             ignoreblobs,
                             "RIVETME");
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,RivetME_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface on top of ME level events.";
}


#endif
