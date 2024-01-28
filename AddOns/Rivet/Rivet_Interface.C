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
#if defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
    HepMC3::GenEvent m_lastevent;
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
#if ! ( defined(USING__MPI) && defined(USING__Rivet_MPI_Merge) )
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

  HepMC3::GenEvent event;
  if (m_usehepmcshort)  m_hepmc2.Sherpa2ShortHepMC(bl, event);
  else                  m_hepmc2.Sherpa2HepMC(bl, event);
  std::vector<HEPMCNS::GenEvent*> subevents(m_hepmc2.GenSubEventList());
  m_hepmc2.AddCrossSection(event, p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
#if defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
  if (m_lastevent.vertices().empty()) {
    m_lastevent=event;
    for (size_t i(0);i<m_lastevent.weights().size();++i) m_lastevent.weights()[i]=0;
  }
  m_hepmc2.AddCrossSection(m_lastevent, p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
  PRINT_VAR(m_lastevent.cross_section());
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
#if defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
  // to be done properly
  const double nomxsec = p_eventhandler->TotalXSMPI().Nominal();
  const double nomxerr = p_eventhandler->TotalErrMPI().Nominal();
#else
  // to be done properly
  const double nomxsec = p_eventhandler->TotalXS().Nominal();
  const double nomxerr = p_eventhandler->TotalErr().Nominal();
#endif

#if defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
    // synchronize analyses among MPI processes
    std::string mynames;
    for (auto& it : m_rivet) {
      std::string out;
      if (it.first.first!="") out+="."+it.first.first;
      if (it.first.second!=0) out+=".j"+ToString(it.first.second);
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
	if (m_rivet.find(key)==m_rivet.end()) {
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
      it.second->collapseEventGroup();
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
#if defined(USING__MPI) && defined(USING__Rivet_MPI_Merge)
      std::vector<double> data;
      it.second->serialize(data);
      mpi->Allreduce(&data[0],data.size(),MPI_DOUBLE,MPI_SUM);
      size_t i(0);
      it.second->deserialize(data,i);
      if (i!=data.size()) THROW(fatal_error,"serialization error");
      if (mpi->Rank()==0)
#endif
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
// end of Rivet_Interface
