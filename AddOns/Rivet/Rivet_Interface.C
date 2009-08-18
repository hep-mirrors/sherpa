#include "SHERPA/Tools/Analysis_Interface.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
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
  typedef std::map<int, AnalysisHandler*> RivetMap;
private:

  std::string m_inpath, m_infile, m_outpath;
  std::vector<std::string> m_analyses;

  size_t m_nevt;
  double m_sum_of_weights;
  bool   m_finished;
  bool   m_splitjetconts;
  
  RivetMap m_rivet;
  HepMC2_Interface m_hepmc2;
  std::vector<btp::code> m_ignoreblobs;

public:

  inline Rivet_Interface(const std::string &inpath,
                         const std::string &infile,
                         const std::string &outpath,
                         const std::vector<btp::code> &ignoreblobs) :
    m_inpath(inpath), m_infile(infile), m_outpath(outpath),
    m_nevt(0), m_sum_of_weights(0.0), m_finished(false),
    m_ignoreblobs(ignoreblobs)
  {
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
      reader.SetFileBegin("BEGIN_RIVET");
      reader.SetFileEnd("END_RIVET");
      
//      m_splitjetconts=reader.GetValue<int>("JETCONTS", 0);
      Log::setLevel("Rivet", reader.GetValue<int>("-l", 20));
      reader.VectorFromFile(m_analyses,"-a");
      m_sum_of_weights=0.0;
      
      for (size_t i=0; i<m_ignoreblobs.size(); ++i) {
        m_hepmc2.Ignore(m_ignoreblobs[i]);
      }
    }
    return true;
  }
  
  AnalysisHandler* GetRivet(int jetcont) {
    RivetMap::iterator it=m_rivet.find(jetcont);
    if (it!=m_rivet.end()) {
      return it->second;
    }
    else {
      std::string out=(jetcont>0 ? m_outpath+".j"+ToString(jetcont) : m_outpath);
      AnalysisHandler* rivet(new AnalysisHandler(out, "", AIDAML));
      rivet->addAnalyses(m_analyses);
      rivet->init();
      m_rivet.insert(std::make_pair(jetcont, rivet));
      return rivet;
    }
  }
  
  
  bool Run(ATOOLS::Blob_List *const bl)
  {
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    double weight((*sp)["Weight"]->Get<double>());
    HepMC::GenEvent event;
    m_hepmc2.Sherpa2HepMC(bl, event, weight);
#ifdef HEPMC_HAS_CROSS_SECTION
    HepMC::GenCrossSection xs;
    xs.set_cross_section(p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
    event.set_cross_section(xs);
#endif
    
    GetRivet(0)->analyze(event);
//    if (m_splitjetconts) {
//      GetRivet(sp->NOutP())->analyze(event);
//    }
    
    ++m_nevt;
    m_sum_of_weights+=weight;
    return true;
  }
  
  
  bool Finish()
  {
    for (RivetMap::iterator it=m_rivet.begin(); it!=m_rivet.end(); ++it) {
//#ifdef USING__RIVET__SETSOW
//      it->second->setSumOfWeights(m_sum_of_weights);
//#endif
      it->second->finalize();
      it->second->commitData();
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
    (args.m_inpath,args.m_infile,outpath, std::vector<btp::code>());
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
    (args.m_inpath,args.m_infile,outpath+".SL", ignoreblobs);
}

void RivetShower_Interface_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface on top of shower level events.";
}

#endif
