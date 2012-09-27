#ifndef EXTRA_XS_Main_Simple_XS_H
#define EXTRA_XS_Main_Simple_XS_H

#include "EXTRA_XS/Main/Process_Group.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace ATOOLS { class Data_Reader;  }
namespace MODEL  { class Model_Base;   }
namespace PDF    { class Remnant_Base; }

namespace EXTRAXS {

  class Cluster_Algorithm;

  class Simple_XS: public Process_Group, public PHASIC::ME_Generator_Base {
  private :

    Cluster_Algorithm *p_cluster;

    std::string m_path, m_file;

    void DrawLogo(std::ostream &ostr);

  public :

    // constructor
    Simple_XS();

    // destructor
    ~Simple_XS();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);
    Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add);
    bool PerformTests();

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs);

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (Process_Base *const proc,const size_t &mode,const double &kt2);

  }; // end of class Simple_XS

} // end of namespace EXTRAXS

#endif

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "EXTRA_XS/Cluster/Cluster_Algorithm.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "EXTRA_XS/Main/Single_Process.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

void Simple_XS::DrawLogo(std::ostream &ostr)
{
}

Simple_XS::Simple_XS(): 
  ME_Generator_Base("Internal"), p_cluster(NULL) 
{
  DrawLogo(std::cout);
}

Simple_XS::~Simple_XS() 
{
  if (p_cluster) delete p_cluster;
}

bool Simple_XS::Initialize(const string &path,const string &file,
			   Model_Base *const model,
			   BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isrhandler)
{
  m_path=path;
  m_file=file;
  p_int->SetBeam(beam); 
  p_int->SetISR(isrhandler);
  return true;
}

Process_Base *Simple_XS::InitializeProcess(const Process_Info &pi, bool add)
{
  size_t n(pi.m_ii.NExternal()+pi.m_fi.NExternal());
  bool oneisgroup(pi.m_ii.IsGroup()||pi.m_fi.IsGroup());
  if (oneisgroup) {
    Process_Group* newxs = new Process_Group();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    if (!newxs->ConstructProcesses()) {
      msg_Debugging()<<METHOD<<"(): Construct failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (add) Add(newxs);
    newxs->SetGenerator(this);
    DEBUG_INFO("Initialized '"<<newxs->Name());
    return newxs;
  }
  else {
    Single_Process* newxs = new Single_Process();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    if (!newxs->Initialize()) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (add) Add(newxs);
    newxs->SetGenerator(this);
    DEBUG_INFO("Initialized '"<<newxs->Name());
    return newxs;
  }
}

bool Simple_XS::PerformTests()
{
  return true;
}
  
void Simple_XS::SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs)
{
  if (p_cluster==NULL) p_cluster = new Cluster_Algorithm();
  p_cluster->SetClusterDefinitions(defs);
}

Cluster_Amplitude *Simple_XS::ClusterConfiguration
(Process_Base *const proc,const size_t &mode,const double &kt2)
{
  if (mode!=2) p_cluster->Cluster(proc->Get<Single_Process>());
  Cluster_Amplitude *ampl(p_cluster->Amplitude());
  ampl->SetMS(this);
  return ampl;
}

namespace PHASIC {

  DECLARE_GETTER(Simple_XS_Getter,"Internal",ME_Generator_Base,ME_Generator_Key);

  ME_Generator_Base *Simple_XS_Getter::operator()(const ME_Generator_Key &key) const
  {
    return new Simple_XS();
  }

  void Simple_XS_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"The internal ME generator"; 
  }

}
