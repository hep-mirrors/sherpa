#include "SHERPA/Tools/Analysis_Interface.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__RIVET

#include "ATOOLS/Math/MathTools.H"
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
#include "SHERPA/Tools/HepMC3_Interface.H"

#include "Rivet/Config/RivetConfig.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenCrossSection.h"


namespace SHERPARIVET {
  typedef std::pair<std::string, int> RivetMapKey;
  typedef std::map<RivetMapKey, Rivet::AnalysisHandler*> Rivet_Map;

  class Rivet_Interface: public SHERPA::Analysis_Interface {
  private:

    std::string m_outpath, m_tag;

    size_t m_nevt, m_histointerval;
    bool   m_finished;
    bool   m_splitjetconts, m_splitSH, m_splitpm, m_splitcoreprocs, m_usehepmcshort;

    int m_loglevel, m_ignorebeams, m_skipweights;
    double m_weightcap, m_nlosmearing;
    std::string m_matchweights, m_unmatchweights, m_nomweight;
    std::vector<std::string> m_analyses;

    Rivet_Map         m_rivet;
    SHERPA::HepMC3_Interface      m_hepmc;
    std::vector<ATOOLS::btp::code> m_ignoreblobs;
    std::map<std::string,size_t>   m_weightidxmap;
#if defined(USING__MPI) && defined(USING__RIVET4)
    HepMC3::GenEvent m_lastevent;
#endif

    Rivet::AnalysisHandler* GetRivet(std::string proc,
                                     int jetcont,
                                     HepMC3::GenEvent* dummyevent=nullptr);
    std::string             GetCoreProc(const std::string& proc);
    std::string             OutputPath(const Rivet_Map::key_type& key);

  public:
    Rivet_Interface(const std::string &outpath,
                    const std::vector<ATOOLS::btp::code> &ignoreblobs,
                    const std::string &tag);
    ~Rivet_Interface();

    bool Init() override;
    bool Run(ATOOLS::Blob_List *const bl) override;
    bool Finish() override;

    void ShowSyntax(const int i) override;
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
  m_splitcoreprocs(false), m_ignoreblobs(ignoreblobs)
{
  // create output path if necessary
  if (m_outpath.rfind('/')!=std::string::npos)
    MakeDir(m_outpath.substr(0,m_outpath.rfind('/')), true);

  // add a MPI rank specific suffix if necessary
#if defined(USING__MPI) && !defined(USING__RIVET4)
  if (mpi->Size()>1)
    m_outpath.insert(m_outpath.length(),"_"+rpa->gen.Variable("RNG_SEED"));
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
                                           int jetcont,
                                           HepMC3::GenEvent* dummyevent)
{
  DEBUG_FUNC(proc<<" "<<jetcont);
  RivetMapKey key = std::make_pair(proc, jetcont);
  Rivet_Map::iterator it=m_rivet.find(key);
  if (it==m_rivet.end()) {
    msg_Debugging()<<"create new "<<key.first<<" "<<key.second<<std::endl;
    m_rivet[key] = new AnalysisHandler();
    m_rivet[key]->addAnalyses(m_analyses);
#ifdef USING__RIVET4
    m_rivet[key]->setCheckBeams(!m_ignorebeams);
    m_rivet[key]->matchWeightNames(m_matchweights);
    m_rivet[key]->unmatchWeightNames(m_unmatchweights);
    m_rivet[key]->setFinalizePeriod(OutputPath(key), m_histointerval);
#else
    m_rivet[key]->setIgnoreBeams(m_ignorebeams);
    m_rivet[key]->selectMultiWeights(m_matchweights);
    m_rivet[key]->deselectMultiWeights(m_unmatchweights);
    m_rivet[key]->setAODump(OutputPath(key), m_histointerval);
#endif
    m_rivet[key]->skipMultiWeights(m_skipweights);
    m_rivet[key]->setNominalWeightName(m_nomweight);
    m_rivet[key]->setWeightCap(m_weightcap);
    m_rivet[key]->setNLOSmearing(m_nlosmearing);
    if (dummyevent)
      m_rivet[key]->init(*dummyevent);
    Log::setLevel("Rivet", m_loglevel);
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
#if !defined(USING__RIVET4)
    if ((m_splitjetconts || m_splitSH || m_splitpm || m_splitcoreprocs)
      && s_variations->HasVariations()) {
      msg_Error()<<"WARNING in "<<METHOD<<":\n"
        <<"   Analysis splitting is combined with on-the-fly variations. Cross\n"
        <<"   sections of variations in split analyses, and hence their\n"
        <<"   normalizations, will not be correct. To fix this, upgrade your\n"
        <<"   Rivet installation to v3.2.0 or later.\n";
    }
#endif

    m_usehepmcshort = s["USE_HEPMC_SHORT"].SetDefault(0).Get<int>();
    if (m_usehepmcshort && m_tag!="RIVET" && m_tag!="RIVETSHOWER") {
      THROW(fatal_error, "Internal error.");
    }

    m_loglevel = s["-l"].SetDefault(1000000).Get<int>();
    m_histointerval = s["HISTO_INTERVAL"].SetSynonyms({"--histo-interval"}).SetDefault(0).Get<size_t>();
    m_ignorebeams = s["IGNORE_BEAMS"].SetSynonyms({"IGNOREBEAMS", "--ignore-beams"}).SetDefault(0).Get<int>();
    m_skipweights = s["SKIP_WEIGHTS"].SetSynonyms({"SKIPWEIGHTS", "--skip-weights"}).SetDefault(0).Get<int>();
    m_weightcap = s["WEIGHT_CAP"].SetSynonyms({"--weight-cap"}).SetDefault(0.0).Get<double>();
    m_nlosmearing = s["NLO_SMEARING"].SetSynonyms({"--nlo-smearing"}).SetDefault(0.0).Get<double>();
    m_matchweights = s["MATCH_WEIGHTS"].SetSynonyms({"--match-weights"}).SetDefault("").Get<std::string>();
    m_unmatchweights = s["UNMATCH_WEIGHTS"].SetSynonyms({"--unmatch-weights"}).SetDefault("").Get<std::string>();
    m_nomweight = s["NOMINAL_WEIGHT"].SetSynonyms({"--nominal-weight"}).SetDefault("").Get<std::string>();
    m_analyses = s["ANALYSES"].SetSynonyms({"ANALYSIS", "-a", "--analyses"})
                              .SetDefault<std::vector<std::string>>({}).GetVector<std::string>();

    // configure HepMC interface
    for (size_t i=0; i<m_ignoreblobs.size(); ++i) {
      m_hepmc.Ignore(m_ignoreblobs[i]);
    }
    m_hepmc.SetHepMCNamedWeights(
        s["USE_HEPMC_NAMED_WEIGHTS"].SetDefault(true).Get<bool>());
    m_hepmc.SetHepMCExtendedWeights(
        s["USE_HEPMC_EXTENDED_WEIGHTS"].SetDefault(false).Get<bool>());
    m_hepmc.SetHepMCTreeLike(
        s["USE_HEPMC_TREE_LIKE"].SetDefault(false).Get<bool>());
  }
  return true;
}

bool Rivet_Interface::Run(ATOOLS::Blob_List *const bl)
{
  DEBUG_FUNC("");

  // get particles and validate momenta
  Particle_List pl=bl->ExtractParticles(part_status::active);
  for (Particle_List::iterator it=pl.begin(); it!=pl.end(); ++it) {
    if ((*it)->Momentum().Nan()) {
      msg_Error()<<METHOD<<" encountered NaN in momentum. Ignoring event:"
                 <<std::endl<<*bl<<std::endl;
      return true;
    }
  }

  // create HepMC (sub)events and add cross section to all of them
  HepMC3::GenEvent event;
  if (m_usehepmcshort)  m_hepmc.Sherpa2ShortHepMC(bl, event);
  else                  m_hepmc.Sherpa2HepMC(bl, event);
  std::vector<HepMC3::GenEvent*> subevents(m_hepmc.GenSubEventList());
  m_hepmc.AddCrossSection(event, p_eventhandler->TotalXS(), p_eventhandler->TotalErr());

#if defined(USING__MPI) && defined(USING__RIVET4)
  if (m_lastevent.vertices().empty()) {
    m_lastevent=event;
    m_lastevent.set_run_info(event.run_info());
    for (size_t i(0);i<m_lastevent.weights().size();++i) m_lastevent.weights()[i]=0;
  }
  m_hepmc.AddCrossSection(m_lastevent, p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
#endif

  // dispatch the events to the main & partial (= split) analysis handlers
  if (subevents.size()) {
    // dispatch subevents to the main analysis handler, then delete them
    for (size_t i(0);i<subevents.size();++i) {
      GetRivet("",0)->analyze(*subevents[i]);
    }
    m_hepmc.DeleteGenSubEventList();
  }
  else {
    // dispatch event to the main analysis handler
    GetRivet("",0)->analyze(event);

    // find the final-state multiplicity used below to bin events into the right
    // partial analysis handlers
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    size_t parts=0;
    if (sp) {
      std::string multi(sp?sp->TypeSpec():"");
      if (multi[3]=='_') multi=multi.substr(2,1);
      else multi=multi.substr(2,2);
      parts=ToType<size_t>(multi);
    }

    // now bin the events into the right partial analysis handlers
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
#if defined(USING__MPI) && defined(USING__RIVET4)
  // synchronize analyses among MPI processes to ensure that all processes
  // have the same analyses set; this is otherwise not guaranteed since we create
  // analyses lazily
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
  mynames.resize(len);
  allnames.resize(len*mpi->Size()+1);
  mpi->Allgather(&mynames[0],len,MPI_CHAR,&allnames[0],len,MPI_CHAR);
  char *catname = new char[len+1];
  for (size_t i(0);i<mpi->Size();++i) {
    snprintf(catname, sizeof(catname),"%s",&allnames[len*i]);
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
      GetRivet(proc,ToType<int>(jets),&m_lastevent);
    }
  }
  delete [] catname;

  // merge Rivet::AnalysisHandlers before finalising
  for (auto& it : m_rivet) {
    std::vector<double> data = it.second->serializeContent(true); //< ensure fixed-length across ranks
    mpi->Reduce(&data[0],data.size(),MPI_DOUBLE,MPI_SUM);
    if (mpi->Rank()==0)
      it.second->deserializeContent(data,(size_t)mpi->Size());
  }
  if (mpi->Rank()==0) {
#endif

#ifdef USING__RIVET4
  GetRivet("",0)->collapseEventGroup();
  // determine weight sums and cross sections when Rivet allows us to properly
  // scale variations in split analyses

  const std::vector<double> sumw = GetRivet("", 0)->weightSumWs();
  const std::vector<std::string> wgtnames = GetRivet("", 0)->weightNames();
#if defined(USING__MPI)
  const auto& xs = p_eventhandler->TotalXSMPI();
  const auto& err = p_eventhandler->TotalErrMPI();
#else
  const auto& xs = p_eventhandler->TotalXS();
  const auto& err = p_eventhandler->TotalErr();
#endif
  std::map<std::string, double> xs_wgts;
  std::map<std::string, double> err_wgts;
  xs.FillVariations(xs_wgts);
  err.FillVariations(err_wgts);
  // At this point, we have a "Nominal" entry (but the Rivet weight name might
  // be different, e.g. an empty string ""), and we might have additional
  // unphysical weights in the Rivet weight sums obtained above. Hence, we make
  // sure that every weight name Rivet reports is filled in our cross section
  // lists. For auxiliary weights, we fill -1.0.
  for (int i {0}; i < wgtnames.size(); i++) {
    if (i == GetRivet("", 0)->defaultWeightIndex()) {
      xs_wgts[wgtnames[i]] = xs.Nominal();
      err_wgts[wgtnames[i]] = err.Nominal();
    }
    else {
      auto it = xs_wgts.find(wgtnames[i]);
      if (it == xs_wgts.end()) {
        xs_wgts[wgtnames[i]] = -1.0;
        err_wgts[wgtnames[i]] = -1.0;
      }
    }
  }

#else
  GetRivet("",0)->finalize();
  // determine the nominal weight sum and the nominal cross section when Rivet
  // does not allow us to properly scale variations in split analyses

  const double nomsumw = GetRivet("",0)->sumW();
#if defined(USING__MPI)
  const double nomxsec = p_eventhandler->TotalXSMPI().Nominal();
  const double nomxerr = p_eventhandler->TotalErrMPI().Nominal();
#else
  const double nomxsec = p_eventhandler->TotalXS().Nominal();
  const double nomxerr = p_eventhandler->TotalErr().Nominal();
#endif

#endif

  // in case additional Rivet instances are used,
  // e.g. for splitting into H+S events/jet multis,
  // these only get to "see" a subset of the events
  // and need to be re-scaled to full cross-section
  const bool needs_rescaling = m_rivet.size() > 1;
  for (auto& it : m_rivet) {
    if (needs_rescaling) {
      // first collapse the event group,
      // then scale the cross-section
      // before finalizing
#ifdef USING__RIVET4
      it.second->collapseEventGroup();
      // determine the weight sums seen by this Rivet run
      const std::vector<double> thissumw = it.second->weightSumWs();
      // calculate and set rescaled cross sections
      std::vector<std::pair<double, double>> this_xs_and_errs (wgtnames.size());
      for (int i {0}; i < wgtnames.size(); i++) {
        if (xs_wgts[wgtnames[i]] == -1.0) {
          // do not rescale unphysical cross sections
          this_xs_and_errs[i] = {-1.0, -1.0};
          continue;
        }
        const double wgtfrac = thissumw[i]/sumw[i];
        this_xs_and_errs[i] = {xs_wgts[wgtnames[i]] * wgtfrac,
          err_wgts[wgtnames[i]] * wgtfrac};
      }
      it.second->setCrossSection(this_xs_and_errs);
#else
      it.second->finalize();
      // determine the weight fraction seen by this Rivet run
      const double wgtfrac = it.second->sumW()/nomsumw;
      // rescale nominal cross-section
      const double thisxs  = nomxsec*wgtfrac;
      const double thiserr = nomxerr*wgtfrac;
      it.second->setCrossSection(thisxs, thiserr);
#endif
    }
    it.second->finalize();

    it.second->writeData(OutputPath(it.first));
  }
#if defined(USING__MPI) && defined(USING__RIVET4)
  } // end of if mpi_rank = 0
#endif
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
  str<<"Rivet interface on top of shower level events";
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
  str<<"Rivet interface on top of ME level events";
}

#endif
// end of Rivet_Interface
