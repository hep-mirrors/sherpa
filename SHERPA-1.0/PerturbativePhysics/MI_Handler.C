#include "MI_Handler.H"

#define USING__Sherpa
#include "Amisic.H"
#include "Data_Read.H"

#ifdef PROFILE__all
#define PROFILE__MI_Handler
#endif
#ifdef PROFILE__MI_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

MI_Handler::MI_Handler(std::string path,std::string file,MODEL::Model_Base *model,
		       BEAM::Beam_Spectra_Handler *beam,PDF::ISR_Handler *isr):
  p_amisic(NULL),
  p_beam(beam),
  p_isr(isr),
  m_type(None),
  m_scalescheme(1),
  m_ycut(1.0e-7)
{
  std::string mihandler;
  ATOOLS::Data_Read *read = new ATOOLS::Data_Read(path+file,true);
  mihandler="None";
  if (read->FileExists()) {
    mihandler=read->GetValue<std::string>("MI_HANDLER",std::string("None"));
    m_scalescheme=read->GetValue<int>("MI_HARD_SCALE",1);
    path+=read->GetValue<std::string>("INPUT_PATH",std::string(""));
    file=read->GetValue<std::string>("INPUT_FILE",file);
  }
  delete read;
  if (mihandler==std::string("Amisic")) {
    p_amisic = new AMISIC::Amisic(model,beam,isr);
    p_amisic->SetInputPath(path);
    p_amisic->SetOutputPath(path);
    p_amisic->SetInputFile(file);
    if (!p_amisic->Initialize()) {
      THROW(fatal_error,"Cannot initialize Amisic.");
    }
    m_ycut=p_amisic->HardBase()->Stop(0);
    m_ycut=ATOOLS::sqr(m_ycut/ATOOLS::rpa.gen.Ecms());
    m_type=Amisic;
  }
}

MI_Handler::~MI_Handler() 
{
  if (p_amisic!=NULL) delete p_amisic;
}

Matrix_Element_Handler *MI_Handler::HardMEHandler()
{
  switch (m_type) {
  case Amisic: return p_amisic->HardMEHandler();
  default    : break;
  }
  return NULL;
}

Matrix_Element_Handler *MI_Handler::SoftMEHandler()
{
  switch (m_type) {
  case Amisic: return p_amisic->SoftMEHandler();
  default    : break;
  }
  return NULL;
}

bool MI_Handler::GenerateHardProcess(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
  switch (m_type) {
  case Amisic: return p_amisic->GenerateHardProcess(blob);
  default    : break;
  }
  return false;
}

bool MI_Handler::GenerateSoftProcess(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
  switch (m_type) {
  case Amisic: return p_amisic->GenerateSoftProcess(blob);
  default    : break;
  }
  return false;
}

bool MI_Handler::GenerateEvent(ATOOLS::Blob_List *bloblist)
{
  switch (m_type) {
  case Amisic: break; // p_amisic->GenerateEvent(bloblist);
  default    : break;
  }
  return false;
}

bool MI_Handler::VetoHardProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
  case Amisic: return p_amisic->VetoHardProcess(blob);
  default    : break;
  }
  return false;
}

void MI_Handler::SetScaleMin(double scalemin,unsigned int i)
{
  switch (m_type) {
  case Amisic:
    p_amisic->HardBase()->SetStop(scalemin,i);
    p_amisic->SoftBase()->SetStart(scalemin,i);
    break;
  default:
    break;
  }
}

void MI_Handler::SetScaleMax(double scalemax,unsigned int i)
{
  switch (m_type) {
  case Amisic:
    p_amisic->HardBase()->SetStart(scalemax,i);
    break;
  default:
    break;
  }
}

double MI_Handler::ScaleMin(unsigned int i)
{
  switch (m_type) {
  case Amisic: return p_amisic->HardBase()->Stop(i);
  default    : break;
  }
  return 0.;
}

double MI_Handler::ScaleMax(unsigned int i)
{
  switch (m_type) {
  case Amisic: return p_amisic->HardBase()->Start(i);
  default    : break;
  }
  return 0.;
}

void MI_Handler::Reset()
{
  switch (m_type) {
  case Amisic: p_amisic->Reset();
  default    : break;
  }
}

void MI_Handler::CleanUp()
{
  switch (m_type) {
  case Amisic: p_amisic->CleanUp();
  default    : break;
  }
}

std::string MI_Handler::MIGenerator() 
{
  return Name();
}

MI_Handler::TypeID MI_Handler::Type() 
{
  return m_type;
}

std::string MI_Handler::Name() 
{
  switch (m_type) {
  case Amisic: return std::string("Amisic");
  case None  : return std::string("None");
  default    : break;
  }
  return std::string("Unknown");
}

unsigned int MI_Handler::NIn()
{
  switch (m_type) {
  case Amisic: return p_amisic->HardXS()->NIn();
  default    : break;
  }
  return 0;
}

unsigned int MI_Handler::NOut()
{
  switch (m_type) {
  case Amisic: return p_amisic->HardXS()->NOut();
  default    : break;
  }
  return 0;
}

PDF::ISR_Handler *MI_Handler::ISRHandler()
{
  return p_isr;
}

