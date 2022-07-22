#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "AMISIC++/Main/Amisic.H"
#include "SHRiMPS/Main/Shrimps.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

MI_Handler::MI_Handler(MODEL::Model_Base *model,
		       PDF::ISR_Handler *isr) :
  p_isr(isr), p_amisic(NULL),
  p_ampl(NULL), p_proc(NULL), p_shower(NULL),
  m_stop(false), m_type(typeID::none), m_name("None")
{
  auto s = Settings::GetMainSettings()["MI_HANDLER"];
  if (!rpa->gen.Beam1().IsHadron() || !rpa->gen.Beam2().IsHadron()) {
    s.OverrideScalar<string>("None");
  }
  m_name = s.SetDefault("Amisic").UseNoneReplacements().Get<string>();
  if (m_name==string("Amisic"))  InitAmisic(model);
  if (m_name==string("Shrimps")) InitShrimps(model);
}

MI_Handler::~MI_Handler() 
{
  if (p_amisic!=NULL) delete p_amisic;
}

void MI_Handler::InitAmisic(MODEL::Model_Base *model)
{
  p_amisic    = new AMISIC::Amisic();
  p_amisic->SetOutputPath(rpa->gen.Variable("SHERPA_RUN_PATH")+"/");
  if (!p_amisic->Initialize(model,p_isr)) {
    msg_Error()<<METHOD<<"(): Cannot initialize MPI generator.\n"
	       <<"   Continue without MPIs and hope for the best.\n";
    delete p_amisic; p_amisic=NULL;
  }
  else m_type = typeID::amisic;
}

void MI_Handler::InitShrimps(MODEL::Model_Base *model)
{
  p_shrimps = new SHRIMPS::Shrimps(p_isr);
  m_type = typeID::shrimps;
}

bool MI_Handler::InitialiseMPIs(const double & scale) 
{
  if (m_type==typeID::amisic) return p_amisic->InitMPIs(scale);
  return true;
}

const Vec4D MI_Handler::SelectPositionForScatter() const {
  if (m_type==typeID::amisic) return p_amisic->SelectPositionForScatter(p_amisic->B());
  return Vec4D(0.,0.,0.,0.);
}

void MI_Handler::SetMaxEnergies(const double & E1,const double & E2) {
  if (m_type==typeID::amisic) p_amisic->SetMaxEnergies(E1,E2); 
}

void MI_Handler::ConnectColours(ATOOLS::Blob * showerblob) {
  if (showerblob) p_remnants->ConnectColours(showerblob);
}

Blob * MI_Handler::GenerateHardProcess()
{
  if (m_type==typeID::amisic) {
    Blob * blob = p_amisic->GenerateScatter();
    if (blob==NULL) m_stop = true;
    return blob;
  }
  return NULL;
}

bool MI_Handler::VetoScatter(Blob *blob)
{
  if (m_type==typeID::amisic) return p_amisic->VetoScatter(blob);
  return true;
}

void MI_Handler::Reset()
{
  m_stop = false;
  if (m_type==typeID::amisic) p_amisic->Reset();
}

void MI_Handler::CleanUp()
{
  m_stop = false;
  if (m_type==typeID::amisic) p_amisic->CleanUp();
}

Cluster_Amplitude * MI_Handler::ClusterConfiguration(Blob * blob)
{
  if (m_type==typeID::amisic)  return p_amisic->ClusterConfiguration(blob);
  if (m_type==typeID::shrimps) return p_shrimps->ClusterConfiguration(blob);
  return NULL;
}

const double MI_Handler::ScaleMin() const
{
  if (m_type==typeID::amisic) return p_amisic->ScaleMin();
  return -1.;
}

const double MI_Handler::ScaleMax() const
{
  if (m_type==typeID::amisic) return p_amisic->ScaleMax();
  return -1.;
}


const double MI_Handler::ImpactParameter() const {
  if (p_amisic) return p_amisic->B();
  return 0.;
}

const bool MI_Handler::IsMinBias() const
{
  if (m_type==typeID::amisic) return p_amisic->IsMinBias();
  return false;
}

void MI_Handler::SetMassMode(const int & massmode) {
  if (m_type==typeID::amisic) p_amisic->SetMassMode(massmode);
}

int MI_Handler::ShiftMasses(Cluster_Amplitude * ampl) {
  if (m_type==amisic) return p_amisic->ShiftMasses(ampl);
  return 0;
}


