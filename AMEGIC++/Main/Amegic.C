#ifndef AMEGIC_Main_Amegic_H
#define AMEGIC_Main_Amegic_H

#include "AMEGIC++/Main/Process_Group.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace AMEGIC {

  class Cluster_Algorithm;

  class Amegic: public Process_Group,
		public PHASIC::ME_Generator_Base {
  private :

    std::string  m_path, m_file;

    MODEL::Model_Base *p_model;

    Cluster_Algorithm *p_cluster;

    std::vector<PHASIC::Process_Base*> m_rsprocs;

    void DrawLogo(std::ostream &ostr);

  public :

    // constructor
    Amegic();

    // destructor
    ~Amegic();

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

  };// end of class Amegic

}// end of namespace AMEGIC

#endif

#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "AMEGIC++/Cluster/Cluster_Algorithm.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;

void Amegic::DrawLogo(std::ostream &ostr)
{
  ostr<<"+-----------------------------------------+\n";
  ostr<<"|   X   X   X XXXX  XXX  XXX  XXX         |\n";
  ostr<<"|  X X  XX XX X    X      X  X     X   X  |\n";
  ostr<<"| X   X X X X XXX  X XXX  X  X    XXX XXX |\n";
  ostr<<"| XXXXX X   X X    X   X  X  X     X   X  |\n";
  ostr<<"| X   X X   X XXXX  XXX  XXX  XXX         |\n";
  ostr<<"+-----------------------------------------+\n";
  ostr<<"| please cite: JHEP 0202:044,2002         |\n";
  ostr<<"+-----------------------------------------+\n";
  rpa.gen.AddCitation
    (1,"Amegic is published under \\cite{Krauss:2001iv}.");
}

Amegic::Amegic(): 
  ME_Generator_Base("Amegic"), p_model(NULL), p_cluster(NULL)
{
  DrawLogo(msg->Info());
  p_testmoms=NULL;
}

Amegic::~Amegic() 
{
  if (p_cluster) delete p_cluster;
}
 
bool Amegic::Initialize(const std::string &path,const std::string &file,
			MODEL::Model_Base *const model,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler)
{
  p_model=model;
  m_path=path;
  m_file=file;
  p_int->SetBeam(beamhandler);
  p_int->SetISR(isrhandler);
  return true;
}

PHASIC::Process_Base *Amegic::InitializeProcess(const PHASIC::Process_Info &pi,
                                                bool add)
{
  if (pi.m_fi.m_nloewtype!=PHASIC::nlo_type::lo) return NULL;
  PHASIC::Process_Base *newxs(NULL);
  size_t nis(pi.m_ii.NExternal()), nfs(pi.m_fi.NExternal());
  size_t nt(pi.m_ii.NTotalExternal()+pi.m_fi.NTotalExternal());
  std::string name(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  Topology top(nis+nfs);
  if (nt>nis+nfs) {
    newxs = new AMEGIC::Process_Group();
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    if (!newxs->Get<AMEGIC::Process_Group>()->
	InitAmplitude(p_model,&top)) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (!newxs->Get<AMEGIC::Process_Group>()->ConstructProcesses()) {
      msg_Debugging()<<METHOD<<"(): Construct failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    newxs->Get<AMEGIC::Process_Group>()->WriteMappingFile();
    msg_Tracking()<<"Initialized '"<<newxs->Name()<<"'\n";
    if (msg_LevelIsTracking()) newxs->Get<AMEGIC::Process_Group>()->PrintProcessSummary();
  }
  else {
    newxs = GetProcess(pi);
    if (!newxs) return NULL;
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    p_testmoms = new Vec4D[newxs->NIn()+newxs->NOut()];
    if (!p_pinfo) {
      p_pinfo = Translate(pi);
      m_nin = newxs->NIn();
      m_flavs.clear();
      for (size_t i=0;i<m_nin;i++) 
	m_flavs.push_back(newxs->Flavours()[i]);
    }
    Phase_Space_Handler::TestPoint(p_testmoms,&newxs->Info());
    newxs->Get<AMEGIC::Process_Base>()->SetTestMoms(p_testmoms);
    newxs->Get<AMEGIC::Process_Base>()->SetPrintGraphs(pi.m_gpath!="");
    if (!newxs->Get<AMEGIC::Process_Base>()->
	InitAmplitude(p_model,&top,m_umprocs,m_errprocs)) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
  }
  if (add) Add(newxs);
  else m_rsprocs.push_back(newxs);
  newxs->SetGenerator(this);
  return newxs;
}

bool Amegic::PerformTests()
{
  bool tests(Process_Group::PerformTests());
  if (NewLibs()) THROW(normal_exit,"New libraries created. Please compile.");
  for (size_t i(0);i<m_rsprocs.size();++i) 
    if (m_rsprocs[i]->Get<AMEGIC::Amegic_Base>()->NewLibs())
      THROW(normal_exit,"New libraries created. Please compile.");
  Minimize();
  return tests;
}

void Amegic::SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs)
{
  if (p_cluster==NULL) p_cluster = new Cluster_Algorithm(this);
  p_cluster->SetClusterDefinitions(defs);
}

Cluster_Amplitude *Amegic::ClusterConfiguration
(PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2)
{
  p_cluster->Cluster(proc->Get<AMEGIC::Process_Base>(),mode,kt2);
  return p_cluster->Amplitude();
}

namespace PHASIC {

  DECLARE_GETTER(Amegic_Getter,"Amegic",ME_Generator_Base,ME_Generator_Key);

  ME_Generator_Base *Amegic_Getter::operator()(const ME_Generator_Key &key) const
  {
    return new Amegic();
  }

  void Amegic_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"The AMEGIC++ ME generator"; 
  }

}

