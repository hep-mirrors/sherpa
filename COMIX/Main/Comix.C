#ifndef COMIX_Main_Comix_H
#define COMIX_Main_Comix_H

#include "COMIX/Main/Process_Group.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace MODEL  { class Model_Base;   }
namespace PDF    { class Remnant_Base; }

namespace COMIX {

  class Single_Process;
  class Cluster_Algorithm;

  class Comix: public Process_Group, public PHASIC::ME_Generator_Base {
  private :

    Cluster_Algorithm *p_cluster;

    std::vector<Single_Process*>       m_umprocs;
    std::vector<PHASIC::Process_Base*> m_rsprocs;

    std::string m_path, m_file;
    int    m_map;
    time_t m_mets;

    void PrintLogo(std::ostream &s);

    void InitVertices(const SP(Model) &model);

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

    PHASIC::Process_Base *GetProcesses();

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
#include "COMIX/Amplitude/Vertex_Base.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <unistd.h>

using namespace COMIX;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

static const PHASIC::nlo_type::code wrongnlotype
(PHASIC::nlo_type::loop|PHASIC::nlo_type::vsub|PHASIC::nlo_type::rsub);

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
  rpa.gen.AddCitation
    (1,"Comix is published under \\cite{Gleisberg:2008fv}.");
}

void Comix::InitVertices(const SP(Model) &model)
{
  Filler_Getter::Getter_List flist(Filler_Getter::GetGetters());
  for (Filler_Getter::Getter_List::const_iterator git(flist.begin());
       git!=flist.end();++git) (*git)->GetObject(&*model);
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n   Implemented currents:\n\n";
    Current_Getter::PrintGetterInfo(msg_Out(),5);
    msg_Out()<<"\n   Implemented vertices:\n\n";
    Vertex_Getter::PrintGetterInfo(msg_Out(),30);
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
  InitModel(model,m_path+m_file);
  if (p_model!=NULL) {
    PrintLogo(msg->Info());
    InitVertices(p_model);
  }
  p_int->SetBeam(beamhandler); 
  p_int->SetISR(isrhandler);
  // init mapping file
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputPath(m_path);
  read.SetInputFile(m_file);
  m_map=read.GetValue<int>("WRITE_MAPPING_FILE",1);
  if (m_map&2) 
    msg_Out()<<om::bold<<METHOD<<"(): "<<om::red
	     <<"Ignoring timestamps on input files."
	     <<om::reset<<std::endl;
  if (m_map>0) {
    struct stat fst;
    std::string fname(rpa.gen.Variable("MODEL_DATA_FILE"));
    if (fname.find('|')!=std::string::npos)
      fname=fname.substr(0,fname.find('|'));
    My_In_File infile(m_path,fname);
    if (infile.Open()) {
      stat((infile.Path()+infile.File()).c_str(),&fst);
      m_mets=fst.st_mtime;
    }
    else {
      msg_Tracking()<<METHOD<<"(): Failed to get timestamp of '"
		    <<rpa.gen.Variable("MODEL_DATA_FILE")<<"'\n";
      m_mets=0;
    }
    fname=rpa.gen.Variable("RUN_DATA_FILE");
    if (fname.find('|')!=std::string::npos)
      fname=fname.substr(0,fname.find('|'));
    infile=My_In_File(m_path,fname);
    if (infile.Open()) {
      stat((infile.Path()+infile.File()).c_str(),&fst);
      m_mets=ATOOLS::Max(fst.st_mtime,m_mets);
    }
    else {
      msg_Tracking()<<METHOD<<"(): Failed to get timestamp of '"
		    <<(m_path+fname)<<"'\n";
    }
  }
  return true;
}

PHASIC::Process_Base *Comix::
InitializeProcess(const PHASIC::Process_Info &pi, bool add)
{
  if (p_model==NULL) return NULL;
  if ((pi.m_fi.m_nloqcdtype&wrongnlotype) || 
      (pi.m_fi.m_nloewtype&wrongnlotype)) return NULL;
  PHASIC::Process_Base *newxs(NULL);
  size_t nis(pi.m_ii.NExternal()), nfs(pi.m_fi.NExternal());
  size_t nt(pi.m_ii.NTotalExternal()+pi.m_fi.NTotalExternal());
  std::string name(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  std::map<std::string,std::string> pmap;
  std::string mapfile(rpa.gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix");
  MakeDir(mapfile,true);
  mapfile+="/"+name+".map";
  struct stat buffer;
  msg_Debugging()<<"checking for '"<<mapfile<<"' ... "<<std::flush;
  if (stat(mapfile.c_str(),&buffer)==-1) {
    msg_Debugging()<<" no mapping info"<<std::endl;
  }
  else {
    if (buffer.st_mtime>m_mets || m_map&2) {
      std::ifstream map(mapfile.c_str());
      if (map.good()) {
	while (!map.eof()) {
	  std::string src, dest;
	  map>>src>>dest;
	  pmap[src]=dest;
	  msg_Debugging()<<" map '"<<src<<"' onto '"<<dest<<"'\n";
	}
      }
    }
  }
  if (nt>nis+nfs) {
    newxs = new Process_Group();
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    if (!newxs->Get<Process_Group>()->Initialize(&pmap,&m_umprocs)) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    newxs->Get<COMIX::Process_Base>()->SetGPath(pi.m_gpath);
    newxs->Get<COMIX::Process_Base>()->SetPSMC(pi.m_psmc);
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
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    newxs->Get<COMIX::Process_Base>()->SetGPath(pi.m_gpath);
    newxs->Get<COMIX::Process_Base>()->SetPSMC(pi.m_psmc);
    if (!newxs->Get<Single_Process>()->Initialize(&pmap,&m_umprocs)) {
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
  remove(mapfile.c_str());
  if (m_map&1) {
    std::ofstream map(mapfile.c_str());
    if (map.good()) {
      for (std::map<std::string,std::string>::const_iterator
	     mit(pmap.begin());mit!=pmap.end();++mit) {
	msg_Debugging()<<" map '"<<mit->first
		       <<"' onto '"<<mit->second<<"'\n";
	map<<mit->first<<" "<<mit->second<<"\n";
      }
    }
  }
  newxs->SetGenerator(this);
  return newxs;
}

PHASIC::Process_Base *Comix::GetProcesses()
{
  return this;
}

bool Comix::PerformTests()
{
  if (!Tests()) return false;
  for (size_t i=0;i<m_rsprocs.size();++i)
    if (!m_rsprocs[i]->Get<COMIX::Process_Base>()->Tests()) return false;
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
