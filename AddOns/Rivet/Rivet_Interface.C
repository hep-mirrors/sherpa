#include "SHERPA/Tools/Analysis_Interface.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "SHERPA/Tools/HepMC2_Interface.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/MyStrStream.H"

#ifdef USING__RIVET
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Logging.hh"

#ifdef USING__HEPMC2__DEFS
#include "HepMC/HepMCDefs.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace Rivet;

class Rivet_Interface: public Analysis_Interface {
  typedef std::map<std::pair<std::string, int>, AnalysisHandler*> RivetMap;
private:

  std::string m_inpath, m_infile, m_outpath, m_tag;
  std::vector<std::string> m_analyses;

  size_t m_nevt;
  double m_sum_of_weights;
  bool   m_finished;
  bool   m_splitjetconts, m_splitcoreprocs, m_usehepmcshort, m_ignorebeams;
  
  RivetMap m_rivet;
  HepMC2_Interface m_hepmc2;
  std::vector<btp::code> m_ignoreblobs;

public:

  inline Rivet_Interface(const std::string &inpath,
                         const std::string &infile,
                         const std::string &outpath,
                         const std::vector<btp::code> &ignoreblobs,
                         const std::string &tag) :
    Analysis_Interface("Rivet"),
    m_inpath(inpath), m_infile(infile), m_outpath(outpath), m_tag(tag),
    m_nevt(0), m_sum_of_weights(0.0), m_finished(false),
    m_ignoreblobs(ignoreblobs)
  {
    if (m_outpath[m_outpath.size()-1]=='/')
      m_outpath=m_outpath.substr(0,m_outpath.size()-1);
  }
  
  
  ~Rivet_Interface()
  {
    if (!m_finished) Finish();
    for (RivetMap::iterator it=m_rivet.begin(); it!=m_rivet.end(); ++it) {
      delete it->second;
    }
  }

  
  bool Init()
  {
    if (m_nevt==0) {
      Data_Reader reader(" ",";","//");
      reader.AddWordSeparator("\t");
      reader.SetAddCommandLine(false);
      reader.SetInputPath(m_inpath);
      std::string infile(m_infile);
      if (infile.find('|')!=std::string::npos)
        infile=infile.substr(0,infile.find('|'));
      reader.SetInputFile(infile);
      reader.AddComment("#");
      reader.SetFileBegin("BEGIN_"+m_tag);
      reader.SetFileEnd("END_"+m_tag);
      reader.AddFileBegin("BEGIN_"+m_tag+"{");
      reader.AddFileEnd("}END_"+m_tag);
      
      m_splitjetconts=reader.GetValue<int>("JETCONTS", 0);
      m_splitcoreprocs=reader.GetValue<int>("SPLITCOREPROCS", 0);
      m_usehepmcshort=reader.GetValue<int>("USE_HEPMC_SHORT", 0);
      m_ignorebeams=reader.GetValue<int>("IGNOREBEAMS", 0);
      Log::setLevel("Rivet", reader.GetValue<int>("-l", 20));
      reader.SetUseGlobalTags(false);
      reader.VectorFromFile(m_analyses,"-a");
      m_sum_of_weights=0.0;
      
      for (size_t i=0; i<m_ignoreblobs.size(); ++i) {
        m_hepmc2.Ignore(m_ignoreblobs[i]);
      }
    }
    return true;
  }
  
  AnalysisHandler* GetRivet(std::string proc, int jetcont) {
    std::pair<std::string, int> key = std::make_pair(proc, jetcont);
    RivetMap::iterator it=m_rivet.find(key);
    if (it!=m_rivet.end()) {
      return it->second;
    }
    else {
      AnalysisHandler* rivet(new AnalysisHandler());
#ifdef USING__RIVET__SETSOW
      rivet->setIgnoreBeams(m_ignorebeams);
#endif
      rivet->addAnalyses(m_analyses);
      m_rivet.insert(std::make_pair(key, rivet));
      return rivet;
    }
  }
  
  std::string GetCoreProc(const std::string& proc) {
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
        if (fl[fl.length()-1]=='b') {
          fl.erase(fl.length()-1,1);
          bar=true;
        }
        else if ((fl[0]=='W' || fl[0]=='H')) {
          if (fl[fl.length()-1]=='-') {
            fl[fl.length()-1]='+';
            bar=true;
          }
        }
        else if (fl[fl.length()-1]=='+') {
          fl[fl.length()-1]='-';
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

    std::vector<Flavour> noresflavs;
    for (size_t i=0; i<nojetflavs.size(); ++i) {
      if (!Flavour(kf_resummed).Includes(nojetflavs[i])) noresflavs.push_back(nojetflavs[i]);
    }

    std::vector<Flavour> finalflavs;
    // start with initial state
    for (size_t i=0; i<2; ++i) {
      if (Flavour(kf_jet).Includes(flavs[i]))
        finalflavs.push_back(Flavour(kf_jet));
      else if (Flavour(kf_resummed).Includes(flavs[i]))
        finalflavs.push_back(Flavour(kf_resummed));
      else
        finalflavs.push_back(flavs[i]);
    }
    // add all non-jet and non-resummed particles
    for (size_t i=0; i<noresflavs.size(); ++i) {
      finalflavs.push_back(noresflavs[i]);
    }
    // add all resummed particles
    for (size_t i=0; i<nojetflavs.size()-noresflavs.size(); ++i) {
      if (finalflavs.size()>3) break;
      finalflavs.push_back(Flavour(kf_resummed));
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
    while (ret[ret.length()-1]=='_') {
      ret.erase(ret.length()-1, 1);
    }

    DEBUG_VAR(ret);
    return ret;
  }
  
  
  bool Run(ATOOLS::Blob_List *const bl)
  {
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    double weight((*sp)["Weight"]->Get<double>());
    HepMC::GenEvent event;
    if (m_usehepmcshort)  m_hepmc2.Sherpa2ShortHepMC(bl, event, weight);
    else                  m_hepmc2.Sherpa2HepMC(bl, event, weight);
    std::vector<HepMC::GenEvent*> subevents(m_hepmc2.GenSubEventList());
#ifdef HEPMC_HAS_CROSS_SECTION
    HepMC::GenCrossSection xs;
    xs.set_cross_section(p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
    event.set_cross_section(xs);
    for (size_t i(0);i<subevents.size();++i) {
      subevents[i]->set_cross_section(xs);
    }
#endif

    if (subevents.size()) {
      for (size_t i(0);i<subevents.size();++i) {
        GetRivet("", 0)->analyze(*subevents[i]);
      }
      m_hepmc2.DeleteGenSubEventList();
    }
    else {
      GetRivet("", 0)->analyze(event);
      if (m_splitjetconts) {
        // temporary fix for powheg & menlops
        // - only works if the core proc is powheg and the others ar LO
        // - nlo w/ powheg emission goes into n_core-analysis
        // - nlo w/o powheg emission goes into (n_core-1)-analysis
        // => still separating these bits
        size_t parts(sp->NOutP());
        if (sp->TypeSpec().find("BVIRS")!=std::string::npos) --parts;
        GetRivet("", parts)->analyze(event);
      }
      if (m_splitcoreprocs) {
        GetRivet(GetCoreProc(sp->TypeSpec()), 0)->analyze(event);
        if (m_splitjetconts) {
          GetRivet(GetCoreProc(sp->TypeSpec()), sp->NOutP())->analyze(event);
        }
      }
    }

    ++m_nevt;
    m_sum_of_weights+=weight;
    return true;
  }
  
  
  bool Finish()
  {
    for (RivetMap::iterator it=m_rivet.begin(); it!=m_rivet.end(); ++it) {
#ifdef USING__RIVET__SETSOW
      it->second->setSumOfWeights(m_sum_of_weights);
#endif
      std::string out=m_outpath;
      if (it->first.first!="") out+="."+it->first.first;
      if (it->first.second!=0) out+=".j"+ToString(it->first.second);
      it->second->finalize();
      it->second->writeData(out+".aida");
    }
    m_finished=true;
    return true;
  }

  
  void ShowSyntax(const int i)
  {
    if (!msg_LevelIsInfo() || i==0) return;
    msg_Out()<<METHOD<<"(): {\n\n"
        <<"   BEGIN_RIVET {\n\n"
        <<"     -a <ana_1> <ana_2>   analyses to run\n";
    msg_Out()<<"\n   } END_RIVET\n\n"
        <<"}"<<std::endl;
  }
  
};// end of class Rivet_Interface


DECLARE_GETTER(Rivet_Interface_Getter,"Rivet",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *Rivet_Interface_Getter::operator()
(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath[outpath.length()-1]=='/') {
    outpath.erase(outpath.length()-1, 1);
  }
  return new Rivet_Interface
    (args.m_inpath,args.m_infile,outpath, std::vector<btp::code>(), "RIVET");
}

void Rivet_Interface_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface";
}


DECLARE_GETTER(RivetShower_Interface_Getter,"RivetShower",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *RivetShower_Interface_Getter::operator()
(const Analysis_Arguments &args) const
{
  std::string outpath=args.m_outpath;
  if (outpath[outpath.length()-1]=='/') {
    outpath.erase(outpath.length()-1, 1);
  }
  std::vector<btp::code> ignoreblobs;
  ignoreblobs.push_back(btp::Fragmentation);
  ignoreblobs.push_back(btp::Hadron_Decay);
  ignoreblobs.push_back(btp::Hadron_Mixing);
  return new Rivet_Interface
    (args.m_inpath,args.m_infile,outpath+".SL", ignoreblobs, "RIVETSHOWER");
}

void RivetShower_Interface_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface on top of shower level events.";
}

#endif
