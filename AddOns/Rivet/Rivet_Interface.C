#include "SHERPA/Tools/Analysis_Interface.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "SHERPA/Tools/HepMC2_Interface.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"

#ifdef USING__RIVET
#include "Rivet/AnalysisHandler.hh"

#ifdef USING__HEPMC2__DEFS
#include "HepMC/HepMCDefs.h"
#ifdef HEPMC_HAS_CROSS_SECTION


using namespace SHERPA;
using namespace ATOOLS;

class Rivet_Interface: public Analysis_Interface {
private:

  std::string m_inpath, m_infile, m_outpath;

  size_t m_nevt;
  bool   m_finished;
  
  Rivet::AnalysisHandler* p_rivet;
  HepMC2_Interface m_hepmc2;

public:

  inline Rivet_Interface(const std::string &inpath,
                         const std::string &infile,
                         const std::string &outpath) :
    m_inpath(inpath), m_infile(infile), m_outpath(outpath),
    m_nevt(0), m_finished(false), p_rivet(NULL)
  {
  }
  
  
  ~Rivet_Interface()
  {
    if (!m_finished) Finish();
    if (p_rivet) delete p_rivet;
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
      std::string fname=reader.GetValue<std::string>("-H", "Rivet");
      
      p_rivet = new Rivet::AnalysisHandler(fname, "", Rivet::AIDAML);
      std::vector<std::string> analyses;
      if (reader.VectorFromFile(analyses,"-a")) {
        p_rivet->addAnalyses(analyses);
      }
      p_rivet->init();
    }
    return true;
  }
  
  
  bool Run(ATOOLS::Blob_List *const bl)
  {
    Blob *sp(bl->FindFirst(btp::Signal_Process));
    double weight((*sp)["Weight"]->Get<double>());
    HepMC::GenEvent event;
    m_hepmc2.Sherpa2HepMC(bl, event, weight);
    HepMC::GenCrossSection xs;
    xs.set_cross_section(p_eventhandler->TotalXS(), p_eventhandler->TotalErr());
    event.set_cross_section(xs);
 
    p_rivet->analyze(event);
    ++m_nevt;
    
    return true;
  }
  
  
  bool Finish()
  {
    p_rivet->finalize();
    p_rivet->commitData();
    m_finished=true;
    return true;
  }

  
  void ShowSyntax(const int i)
  {
    if (!msg_LevelIsInfo() || i==0) return;
    msg_Out()<<METHOD<<"(): {\n\n"
        <<"   BEGIN_RIVET {\n\n"
        <<"     -H <filename>        output file name\n"
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
  return new Rivet_Interface
    (args.m_inpath,args.m_infile,args.m_outpath);
}

void Rivet_Interface_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"Rivet interface";
}

#endif
#endif
#endif
