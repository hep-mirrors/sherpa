#include "MI_Handler.H"

#include "Data_Read.H"
#include "Lund_Interface.H"
#include "Lund_Wrapper.H"

using namespace SHERPA;

MI_Handler::MI_Handler(std::string path,std::string file,MODEL::Model_Base *model,
		       BEAM::Beam_Spectra_Handler *beam,PDF::ISR_Handler *isr):
  p_amisic(NULL),
  p_beam(beam),
  p_isr(isr),
  m_type(None)
{
  std::string mihandler, beamfile;
  ATOOLS::Data_Read *read = new ATOOLS::Data_Read(path+file);
  mihandler=read->GetValue<std::string>("MI_HANDLER",std::string("Amisic"));
  path+=read->GetValue<std::string>("INPUT_PATH",std::string(""));
  file=read->GetValue<std::string>("INPUT_FILE",file);
  beamfile=read->GetValue<std::string>("BEAM_DATA_FILE","Beam.dat");
  delete read;
  if (mihandler==std::string("Amisic")) {
    p_amisic = new AMISIC::Amisic(model,beam,isr);
    p_amisic->SetInputPath(path);
    p_amisic->SetOutputPath(path);
    p_amisic->SetInputFile(file);
    p_amisic->Initialize();
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
  case Amisic:
    return p_amisic->HardMEHandler();
    break;
  default:
    break;
  }
  return false;
}

Matrix_Element_Handler *MI_Handler::SoftMEHandler()
{
  switch (m_type) {
  case Amisic:
    return p_amisic->SoftMEHandler();
    break;
  default:
    break;
  }
  return false;
}

bool MI_Handler::GenerateHardProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
  case Amisic:
    return p_amisic->GenerateHardProcess(blob);
    break;
  default:
    break;
  }
  return false;
}

bool MI_Handler::GenerateSoftProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
  case Amisic:
    return p_amisic->GenerateSoftProcess(blob);
    break;
  default:
    break;
  }
  return false;
}

void MI_Handler::SameHardProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
  case Amisic:
    p_amisic->SameHardProcess(blob);
    break;
  default:
    break;
  }
}

void MI_Handler::SameSoftProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
  case Amisic:
    p_amisic->SameSoftProcess(blob);
    break;
  default:
    break;
  }
}

void MI_Handler::SetScaleMin(double scalemin,unsigned int i)
{
  switch (m_type) {
  case Amisic:
    p_amisic->HardBase()->SetStop(scalemin,i);
//    p_amisic->SoftBase()->SetStop(scalemin,i);
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

void MI_Handler::Reset()
{
  switch (m_type) {
  case Amisic:
    p_amisic->Reset();
    break;
  default:
    break;
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
  case Amisic:
    return p_amisic->HardXS()->Nin();
    break;
  default:
    break;
  }
  return 0;
}

unsigned int MI_Handler::NOut()
{
  switch (m_type) {
  case Amisic:
    return 2;
    break;
  default:
    break;
  }
  return 0;
}

PDF::ISR_Handler *MI_Handler::ISRHandler()
{
  return p_isr;
}

