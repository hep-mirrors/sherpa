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
		       PDF::ISR_Handler *isr, YFS::YFS_Handler *yfs,
                       REMNANTS::Remnant_Handler * remnant_handler) :
  p_isr(isr), p_yfs(yfs), m_id(p_isr->Id()), p_remnants(remnant_handler),
  p_amisic(NULL), p_shrimps(NULL), p_ampl(NULL),
  p_proc(NULL), p_shower(NULL), m_on(true), m_stop(false),
  m_firstrescatter((m_id==PDF::isr::bunch_rescatter) ? true : false),
  m_type(typeID::none), m_name("None")
{
  Settings& s = Settings::GetMainSettings();
  m_name      = s["MI_HANDLER"].SetDefault("Amisic").UseNoneReplacements().Get<string>();
  string scm  = s["SOFT_COLLISIONS"].SetDefault("None").UseNoneReplacements().Get<string>();
  if (m_id==PDF::isr::bunch_rescatter) {
    string resc = s["BEAM_RESCATTERING"].Get<string>();
    scm = m_name = resc;
  }
  if (isr->Mode() != PDF::isrmode::hadron_hadron || m_name=="None") {
    m_name = "None";
    m_on   = false;
  }
  else {
    if (m_name==string("Amisic"))  InitAmisic(model);
    if ((scm==string("Shrimps") && p_amisic==NULL) ||
	m_name==string("Shrimps")) InitShrimps(model);
  }
  msg_Out()<<"Multiple interactions initialized\n"
           <<"  id = "<<m_id<<", name = "<<m_name<<", type = "<<m_type<<"\n";
}

MI_Handler::~MI_Handler() 
{
  if (p_amisic!=NULL)  { delete p_amisic;  p_amisic  = NULL; }
  if (p_shrimps!=NULL) { delete p_shrimps; p_shrimps = NULL; }
}

void MI_Handler::InitAmisic(MODEL::Model_Base *model)
{
  p_amisic    = new AMISIC::Amisic();
  p_amisic->SetOutputPath(rpa->gen.Variable("SHERPA_RUN_PATH")+"/");
  if (!p_amisic->Initialize(model,p_isr,p_yfs,p_remnants)) {
    msg_Error()<<METHOD<<"(): Cannot initialize MPI generator. \n"
	       <<"Continue without MPIs and hope for the best.\n";
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
  if (m_type==typeID::amisic) return p_amisic->InitMPIs(p_isr, scale);
  return true;
}

const Vec4D MI_Handler::SelectPositionForScatter() const {
  if (m_type==typeID::amisic) return p_amisic->SelectPositionForScatter(p_amisic->B());
  return Vec4D(0.,0.,0.,0.);
}

void MI_Handler::SetMaxEnergies(const double & E1,const double & E2) {
  if (m_type==typeID::amisic)  p_amisic->SetMaxEnergies(E1,E2);
  if (m_type==typeID::shrimps) p_shrimps->SetMaxEnergies(E1,E2);
}

void MI_Handler::ConnectColours(ATOOLS::Blob * showerblob) {
  if (!m_firstrescatter && showerblob) p_remnants->ConnectColours(showerblob);
}

Blob * MI_Handler::GenerateHardProcess()
{
  Blob * blob = NULL;
  if (m_type==typeID::amisic)  blob = p_amisic->GenerateScatter();
  if (m_type==typeID::shrimps) blob = p_shrimps->GenerateEvent();
  if (blob==NULL) m_stop = true;
  m_firstrescatter = false;
  return blob;
}

bool MI_Handler::VetoScatter(Blob *blob)
{
  // Method not yet implemented in either Shrimps or Amisic
  return false;
}

void MI_Handler::Reset()
{
  m_stop = false;
  if (m_type==typeID::amisic) p_amisic->Reset();
  for (short unsigned int i=0;i<2;++i) {
    p_remnants->GetRemnant(i)->Reset(m_id);
    p_isr->ResetRescaleFactor(i);
    p_isr->Reset(i);
  }
}

void MI_Handler::CleanUp()
{
  m_stop           = false;
  m_firstrescatter = (m_id==PDF::isr::bunch_rescatter) ? true : false;
  if (m_type==typeID::amisic)  p_amisic->CleanUp();
  if (m_type==typeID::shrimps) p_shrimps->CleanUp();
}

Cluster_Amplitude * MI_Handler::ClusterConfiguration(Blob * blob)
{
  if (m_type==typeID::amisic)  return p_amisic->ClusterConfiguration(blob);
  if (m_type==typeID::shrimps) return p_shrimps->ClusterConfiguration(blob);
  return NULL;
}

const double MI_Handler::ScaleMin() const
{
  if (m_type==typeID::amisic)  return p_amisic->ScaleMin();
  if (m_type==typeID::shrimps) return p_shrimps->ScaleMin();
  return -1.;
}

const double MI_Handler::ScaleMax() const
{
  if (m_type==typeID::amisic) return p_amisic->ScaleMax();
  if (m_type==typeID::shrimps) return p_shrimps->ScaleMax();
  return -1.;
}


const double MI_Handler::ImpactParameter() const {
  if (m_type==typeID::amisic)  return p_amisic->B();
  if (m_type==typeID::shrimps) return p_shrimps->B();
  return 0.;
}

const bool MI_Handler::IsMinBias() const
{
  if (m_type==typeID::amisic)  return p_amisic->IsMinBias();
  if (m_type==typeID::shrimps) return p_shrimps->IsMinBias();
  return false;
}

void MI_Handler::SetMassMode(const int & massmode) {
  if (m_type==typeID::amisic) p_amisic->SetMassMode(massmode);
}

int MI_Handler::ShiftMasses(Cluster_Amplitude * ampl) {
  if (m_type==amisic) return p_amisic->ShiftMasses(ampl);
  return 0;
}


