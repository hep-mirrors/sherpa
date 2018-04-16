#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "AMISIC++/Main/Amisic.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

MI_Handler::MI_Handler(string path,string file,MODEL::Model_Base *model,
		       PDF::ISR_Handler *isr) :
  p_isr(isr),p_amisic(NULL),p_ampl(NULL),p_proc(NULL),p_shower(NULL),
  m_stop(false),m_type(None),m_name("None")
{
  if (!rpa->gen.Beam1().IsHadron() || !rpa->gen.Beam2().IsHadron()) return;
  Default_Reader read;
  read.SetInputPath(path);
  read.SetInputFile(file);
  string mihandler=read.GetValue<string>("MI_HANDLER","Amisic");
  if (mihandler==string("Amisic")) InitAmisic(path,&read,model);
}

MI_Handler::~MI_Handler() 
{
  if (p_amisic!=NULL) delete p_amisic;
}


void MI_Handler::InitAmisic(string & path,Default_Reader *const dr,
			    MODEL::Model_Base *model)
{
  path       += dr->GetValue<string>("INPUT_PATH","");
  string file = dr->GetValue<string>("INPUT_FILE","");
  p_amisic    = new AMISIC::Amisic();
  p_amisic->SetInputPath(path);
  p_amisic->SetOutputPath(rpa->gen.Variable("SHERPA_RUN_PATH")+"/");
  p_amisic->SetInputFile(file);
  if (!p_amisic->Initialize(dr,model,p_isr)) {
    msg_Error()<<METHOD<<"(): Cannot initialize MPI generator. "
	       <<"Continue without.\n";
    delete p_amisic; p_amisic=NULL;
  }
  m_type=Amisic;
  m_name="AMISIC";
}

bool MI_Handler::InitialiseMPIs(const double & scale) 
{
  if (m_type==Amisic) {
    p_amisic->SetMassMode(1);
    p_amisic->SetMaxScale(scale);
    p_amisic->SetB();
    if (p_amisic->VetoEvent(scale)) {
      m_stop = true;
      return false;
    }
  }
  return true;
}

void MI_Handler::SetMaxEnergies(const double & E1,const double & E2) {
  if (m_type==Amisic) p_amisic->SetMaxEnergies(E1,E2); 
}

void MI_Handler::ConnectColours(ATOOLS::Blob * showerblob) {
  p_remnants->ConnectColours(showerblob);
}

Blob * MI_Handler::GenerateHardProcess()
{
  if (m_type==Amisic) {
    Blob * blob = p_amisic->GenerateScatter();
    if (blob==NULL) m_stop = true;
    return blob;
  }
  return NULL;
}

bool MI_Handler::VetoScatter(Blob *blob)
{
  if (m_type==Amisic) return p_amisic->VetoScatter(blob);
  return true;
}

void MI_Handler::Reset()
{
  m_stop = false;
  if (m_type==Amisic) p_amisic->Reset();
}

void MI_Handler::CleanUp()
{
  m_stop = false;
  if (m_type==Amisic) p_amisic->CleanUp();
}

Cluster_Amplitude * MI_Handler::ClusterConfiguration(Blob * blob)
{
  if (m_type==Amisic) return p_amisic->ClusterConfiguration(blob);
  return NULL;
}

const double MI_Handler::ScaleMin() const
{
  if (m_type==Amisic) return p_amisic->ScaleMin();
  return -1.;
}

const double MI_Handler::ScaleMax() const
{
  if (m_type==Amisic) return p_amisic->ScaleMax();
  return -1.;
}

void MI_Handler::SetMassMode(const int & massmode) {
  if (m_type==Amisic) p_amisic->SetMassMode(massmode);
}

int MI_Handler::ShiftMasses(Cluster_Amplitude * ampl) {
  if (m_type==Amisic) return p_amisic->ShiftMasses(ampl);
  return 0;
}


