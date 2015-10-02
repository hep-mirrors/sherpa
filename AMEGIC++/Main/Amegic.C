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

    MODEL::Model_Base *p_mmodel;
    Amegic_Model      *p_amodel;

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
    int PerformTests();
    bool NewLibraries();

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs);

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const ATOOLS::Vec4D_Vector &p,
     const size_t &mode);

  };// end of class Amegic

}// end of namespace AMEGIC

#endif

#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "AMEGIC++/Cluster/Cluster_Algorithm.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/UFO/UFO_Model.H"

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
  rpa->gen.AddCitation
    (1,"Amegic is published under \\cite{Krauss:2001iv}.");
}

Amegic::Amegic(): 
  ME_Generator_Base("Amegic"), p_mmodel(NULL), p_amodel(NULL), p_cluster(NULL)
{
  DrawLogo(msg->Info());
  p_testmoms=NULL;
  p_gen=this;
}

Amegic::~Amegic() 
{
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  if (p_cluster) delete p_cluster;
  delete p_amodel;
}
 
bool Amegic::Initialize(const std::string &path,const std::string &file,
			MODEL::Model_Base *const model,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler)
{
  Data_Reader read(" ",";","#","=");
  read.AddWordSeparator("\t");
  read.SetInputPath(m_path);
  read.SetInputFile(m_file);
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model) &&
      !read.GetValue<int>("AMEGIC_ALLOW_UFO", 0))
    THROW(fatal_error, "AMEGIC can only be used in built-in models. Please use Comix for UFO models.");
  p_mmodel=model;
  p_amodel = new Amegic_Model(model);
  m_path=path;
  m_file=file;
  p_int->SetBeam(beamhandler);
  p_int->SetISR(isrhandler);
  SetPSMasses(&read);
  double alpha=read.GetValue<double>("AMEGIC_TCHANNEL_ALPHA",0.9);
  rpa->gen.SetVariable("AMEGIC_TCHANNEL_ALPHA",ToString(alpha));
  double salpha=read.GetValue<double>("AMEGIC_SCHANNEL_ALPHA",0.75);
  rpa->gen.SetVariable("AMEGIC_SCHANNEL_ALPHA",ToString(salpha));
  double eps=read.GetValue<double>("AMEGIC_CHANNEL_EPSILON",0.0);
  rpa->gen.SetVariable("AMEGIC_CHANNEL_EPSILON",ToString(eps));
  int gauge(read.GetValue<int>("AMEGIC_DEFAULT_GAUGE",1));
  AMEGIC::Process_Base::SetGauge(gauge);
  if (gauge!=10) msg_Info()<<METHOD<<"(): Set gauge "<<gauge<<"."<<std::endl;

  MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process",true);
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  return true;
}

PHASIC::Process_Base *Amegic::InitializeProcess(const PHASIC::Process_Info &pi,
                                                bool add)
{
  if (pi.m_fi.m_nloewtype!=PHASIC::nlo_type::lo) return NULL;
  PHASIC::Process_Base *newxs(NULL);
  size_t nis(pi.m_ii.NExternal()), nfs(pi.m_fi.NExternal());
  std::string name(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  Topology top(nis+nfs);
  bool oneisgroup(pi.m_ii.IsGroup()||pi.m_fi.IsGroup());
  if (oneisgroup) {
    newxs = new AMEGIC::Process_Group();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    if (!newxs->Get<AMEGIC::Process_Group>()->
	InitAmplitude(p_amodel,&top)) {
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
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    p_testmoms = new Vec4D[newxs->NIn()+newxs->NOut()];
    if (!p_pinfo) {
      p_pinfo = Translate(pi);
      m_nin = newxs->NIn();
      m_flavs.clear();
      for (size_t i=0;i<m_nin;i++) 
	m_flavs.push_back(newxs->Flavours()[i]);
    }
    Phase_Space_Handler::TestPoint(p_testmoms,&newxs->Info(),this);
    Vec4D sum;
    Poincare lab(Vec4D(sqrt(10.0),0.0,0.0,1.0));
    msg_Debugging()<<"After boost:\n";
    for (size_t i(0);i<nis+nfs;++i) {
      lab.Boost(p_testmoms[i]);
      sum+=i<m_nin?-p_testmoms[i]:p_testmoms[i];
      msg_Debugging()<<"  p["<<i<<"] = "<<p_testmoms[i]<<"\n";
    }
    msg_Debugging()<<"} -> sum = "<<sum<<"\n";
    Poincare rot(Vec4D::ZVEC,Vec4D(sqrt(14.0),1.0,2.0,3.0));
    msg_Debugging()<<"After rotation {\n";
    for (size_t i(0);i<nis+nfs;++i) {
      rot.Rotate(p_testmoms[i]);
      sum+=i<m_nin?-p_testmoms[i]:p_testmoms[i];
      msg_Debugging()<<"  p["<<i<<"] = "<<p_testmoms[i]<<"\n";
    }
    msg_Debugging()<<"} -> sum = "<<sum<<"\n";
    newxs->Get<AMEGIC::Process_Base>()->SetTestMoms(p_testmoms);
    newxs->Get<AMEGIC::Process_Base>()->SetPrintGraphs(pi.m_gpath);
    My_In_File::ExecDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/","begin");
    if (!newxs->Get<AMEGIC::Process_Base>()->
	InitAmplitude(p_amodel,&top,m_umprocs,m_errprocs)) {
      My_In_File::ExecDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/","commit");
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    My_In_File::ExecDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/","commit");
  }
  if (add) Add(newxs);
  else m_rsprocs.push_back(newxs);
  newxs->SetGenerator(this);
  return newxs;
}

int Amegic::PerformTests()
{
  int tests(Process_Group::PerformTests());
  if (NewLibs()) THROW(normal_exit,"New libraries created. Please compile.");
  for (size_t i(0);i<m_rsprocs.size();++i) 
    if (m_rsprocs[i]->Get<AMEGIC::Amegic_Base>()->NewLibs())
      THROW(normal_exit,"New libraries created. Please compile.");
  Minimize();
  return tests;
}

bool Amegic::NewLibraries()
{
  if (NewLibs()) return true;
  for (size_t i(0);i<m_rsprocs.size();++i)
    if (m_rsprocs[i]->Get<AMEGIC::Amegic_Base>()->NewLibs()) return true;
  return false;
}

void Amegic::SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs)
{
  if (p_cluster==NULL) p_cluster = new Cluster_Algorithm(this);
  p_cluster->SetClusterDefinitions(defs);
}

Cluster_Amplitude *Amegic::ClusterConfiguration
(PHASIC::Process_Base *const proc,const Vec4D_Vector &p,
 const size_t &mode)
{
  p_cluster->Cluster(proc->Get<AMEGIC::Process_Base>(),mode);
  return p_cluster->Amplitude();
}

DECLARE_GETTER(Amegic,"Amegic",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,Amegic>::
operator()(const ME_Generator_Key &key) const
{
  return new Amegic();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,Amegic>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The AMEGIC++ ME generator"; 
}

