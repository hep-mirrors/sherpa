#ifndef COMIX_Main_Comix_H
#define COMIX_Main_Comix_H

#include "COMIX/Main/Process_Group.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__MPI
#include "mpi.h"
#endif

namespace MODEL  { class Model_Base;   }
namespace PDF    { class Remnant_Base; }

namespace COMIX {

  class Single_Process;
  class Cluster_Algorithm;

  class Comix: public Process_Group, public PHASIC::ME_Generator_Base {
  private :

    Cluster_Algorithm *p_cluster;

    std::vector<std::vector<Single_Process*> > m_umprocs;
    std::vector<PHASIC::Process_Base*>         m_rsprocs;

    std::string m_path, m_file;
    int    m_break;
    time_t m_mets;

    void PrintLogo(std::ostream &s);
    void PrintVertices();

  public :

    // constructor
    Comix();

    // destructor
    ~Comix();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
                    MODEL::Model_Base *const model,
                    BEAM::Beam_Spectra_Handler *const beamhandler,
                    PDF::ISR_Handler *const isrhandler);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi,
                                            bool add);
    bool PerformTests();

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs);

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2);

  }; // end of class Comix

} // end of namespace COMIX

#endif

#include "COMIX/Main/Single_Process.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Cluster/Cluster_Algorithm.H"
#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/CXXFLAGS.H"

using namespace COMIX;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Comix::Comix(): 
  ME_Generator_Base("Comix"), p_cluster(NULL)
{
}

Comix::~Comix() 
{
  if (p_cluster) delete p_cluster;
}

#define RED(ARG) om::red<<ARG<<om::reset
#define GREEN(ARG) om::green<<ARG<<om::reset
#define BLUE(ARG) om::blue<<ARG<<om::reset
#define YELLOW(ARG) om::brown<<ARG<<om::reset
#define BLACK(ARG) ARG

void Comix::PrintLogo(std::ostream &s)
{
  s<<"+----------------------------------+\n";
  s<<"|                                  |\n";
  s<<"|      "<<RED("CCC")<<"  "<<GREEN("OOO")<<"  "
   <<BLUE("M")<<"   "<<BLUE("M")<<" "<<BLACK("I")<<" "
   <<YELLOW("X")<<"   "<<YELLOW("X")<<"     |\n";
  s<<"|     "<<RED("C")<<"    "<<GREEN("O")<<"   "
   <<GREEN("O")<<" "<<BLUE("MM")<<" "<<BLUE("MM")
   <<" "<<BLACK("I")<<"  "<<YELLOW("X")<<" "
   <<YELLOW("X")<<"      |\n";
  s<<"|     "<<RED("C")<<"    "<<GREEN("O")
   <<"   "<<GREEN("O")<<" "<<BLUE("M")<<" "
   <<BLUE("M")<<" "<<BLUE("M")<<" "<<BLACK("I")
   <<"   "<<YELLOW("X")<<"       |\n";
  s<<"|     "<<RED("C")<<"    "<<GREEN("O")
   <<"   "<<GREEN("O")<<" "<<BLUE("M")<<"   "
   <<BLUE("M")<<" "<<BLACK("I")<<"  "
   <<YELLOW("X")<<" "<<YELLOW("X")<<"      |\n";
  s<<"|      "<<RED("CCC")<<"  "<<GREEN("OOO")
   <<"  "<<BLUE("M")<<"   "<<BLUE("M")<<" "
   <<BLACK("I")<<" "<<YELLOW("X")<<"   "
   <<YELLOW("X")<<"     |\n";
  s<<"|                                  |\n";
  s<<"+==================================+\n";
  s<<"|  Color dressed  Matrix Elements  |\n";
  s<<"|     http://comix.freacafe.de     |\n";
  s<<"|   please cite  JHEP12(2008)039   |\n";
  s<<"+----------------------------------+\n";
#ifdef USING__Threading
  s<<"Comix was compiled for multithreading.\n";
#endif
  rpa->gen.AddCitation
    (1,"Comix is published under \\cite{Gleisberg:2008fv}.");
}

void Comix::PrintVertices()
{
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n   Implemented currents:\n\n";
    Current_Getter::PrintGetterInfo(msg_Out(),10);
    msg_Out()<<"\n   Implemented lorentz calculators:\n\n";
    LC_Getter::PrintGetterInfo(msg_Out(),10);
    msg_Out()<<"\n   Implemented color calculators:\n\n";
    CC_Getter::PrintGetterInfo(msg_Out(),10);
    msg_Out()<<"\n}\n";
  }
}

bool Comix::Initialize(const std::string &path,const std::string &file,
		       MODEL::Model_Base *const model,
		       BEAM::Beam_Spectra_Handler *const beamhandler,
		       PDF::ISR_Handler *const isrhandler) 
{
  m_path=path;
  m_file=file;
  p_model=model;
  PrintLogo(msg->Info());
  PrintVertices();
  p_int->SetBeam(beamhandler); 
  p_int->SetISR(isrhandler);
  // init mapping file
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputPath(m_path);
  read.SetInputFile(m_file);
  m_break=read.GetValue<int>("ONLY_MAPPING_FILE",0);
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()==0)
#endif
  MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")
	  +"/Process/Comix",true);
  return true;
}

PHASIC::Process_Base *Comix::
InitializeProcess(const PHASIC::Process_Info &pi, bool add)
{
  if (p_model==NULL) return NULL;
  m_umprocs.push_back(std::vector<Single_Process*>());
  PHASIC::Process_Base *newxs(NULL);
  size_t nis(pi.m_ii.NExternal()), nfs(pi.m_fi.NExternal());
  size_t nt(pi.m_ii.NTotalExternal()+pi.m_fi.NTotalExternal());
  std::map<std::string,std::string> pmap;
  if (nt>nis+nfs) {
    newxs = new Process_Group();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    if (!newxs->Get<Process_Group>()->Initialize(&pmap,&m_umprocs.back())) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    newxs->Get<COMIX::Process_Base>()->SetGPath(pi.m_gpath);
    if (!newxs->Get<PHASIC::Process_Group>()->ConstructProcesses()) {
      msg_Debugging()<<METHOD<<"(): Construct failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    msg_Tracking()<<"Initialized '"<<newxs->Name()<<"'\n";
  }
  else {
    newxs = new Single_Process();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    newxs->Get<COMIX::Process_Base>()->SetGPath(pi.m_gpath);
    if (!newxs->Get<Single_Process>()->Initialize(&pmap,&m_umprocs.back())) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (!newxs->Get<Single_Process>()->MapProcess())
      if (!msg_LevelIsTracking()) msg_Info()<<"."<<std::flush;
  }
  if (add) Add(newxs);
  else m_rsprocs.push_back(newxs);
  return newxs;
}

bool Comix::PerformTests()
{
  if (!Tests()) return false;
  for (size_t i=0;i<m_rsprocs.size();++i)
    if (!m_rsprocs[i]->Get<COMIX::Process_Base>()->Tests()) return false;
  if (m_break) {
    msg_Out()<<std::endl;
    THROW(normal_exit,"Mapping files created. Stop upon request.");
  }
  return true;
}

void Comix::SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs)
{
  if (p_cluster==NULL) p_cluster = new Cluster_Algorithm(this);
  p_cluster->SetClusterDefinitions(defs);
}

Cluster_Amplitude *Comix::ClusterConfiguration
(PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2)
{
  p_cluster->Cluster(proc->Get<COMIX::Single_Process>(),mode,kt2);
  return p_cluster->GetAmplitude();
}

namespace PHASIC {

  DECLARE_GETTER(Comix_Getter,"Comix",ME_Generator_Base,ME_Generator_Key);

  ME_Generator_Base *Comix_Getter::operator()(const ME_Generator_Key &key) const
  {
    return new Comix();
  }

  void Comix_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"The Comix ME generator"; 
  }

}
