#include "MI_Handler.H"
#include "Data_Read.H"

#include "Matrix_Element_Handler.H"

#ifdef USING__Amisic
#include "Amisic.H"
#endif

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
		       BEAM::Beam_Spectra_Handler *beam,PDF::ISR_Handler *isr) :
  p_beam(beam),
  p_isr(isr),
#ifdef USING__Amisic
  p_amisic(NULL),
#endif
  p_hardmehandler(NULL), p_softmehandler(NULL),
  m_type(None),
  m_scalescheme(1),
  m_ycut(1.0e-7)
{
#ifdef USING__Amisic
  std::string mihandler="None";
  ATOOLS::Data_Read read(path+file,true);
  if (read.FileExists()) {
    mihandler=read.GetValue<std::string>("MI_HANDLER",std::string("Amisic"));
    m_scalescheme=read.GetValue<int>("MI_HARD_SCALE",1);
    path+=read.GetValue<std::string>("INPUT_PATH",std::string(""));
    file=read.GetValue<std::string>("INPUT_FILE",file);
  }
  if (!ATOOLS::rpa.gen.Beam1().IsHadron() ||
      !ATOOLS::rpa.gen.Beam2().IsHadron()) mihandler="None";
  if (mihandler==std::string("Amisic")) {
    p_amisic = new AMISIC::Amisic(model,beam,isr);
    p_amisic->SetInputPath(path);
    p_amisic->SetOutputPath(path);
    p_amisic->SetInputFile(file);
    if (!p_amisic->Initialize()) {
      THROW(fatal_error,"Cannot initialize Amisic.");
    }
    p_hardmehandler = new Matrix_Element_Handler();
    p_hardmehandler->SetXS((EXTRAXS::Simple_XS*)p_amisic->HardBase()->XS());
    p_hardmehandler->SetUseSudakovWeight(p_amisic->HardBase()->JetVeto());
    p_softmehandler = new Matrix_Element_Handler();
    p_softmehandler->SetXS((EXTRAXS::Simple_XS*)p_amisic->SoftBase()->XS());
    p_softmehandler->SetUseSudakovWeight(p_amisic->SoftBase()->JetVeto());
    m_ycut=p_amisic->HardBase()->Stop(0);
    m_ycut=ATOOLS::sqr(m_ycut/ATOOLS::rpa.gen.Ecms());
    m_type=Amisic;
  }
#endif
}

MI_Handler::~MI_Handler() 
{
#ifdef USING__Amisic
  if (p_amisic!=NULL) delete p_amisic;
#endif
  if (p_hardmehandler!=NULL) delete p_hardmehandler;
  if (p_softmehandler!=NULL) delete p_softmehandler;
}

bool MI_Handler::GenerateHardProcess(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->GenerateHardProcess(blob);
#endif
  default    : break;
  }
  return false;
}

bool MI_Handler::GenerateSoftProcess(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->GenerateSoftProcess(blob);
#endif
  default    : break;
  }
  return false;
}

bool MI_Handler::GenerateEvent(ATOOLS::Blob_List *bloblist)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: break; // p_amisic->GenerateEvent(bloblist);
#endif
  default    : break;
  }
  return false;
}

bool MI_Handler::VetoHardProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->VetoHardProcess(blob);
#endif
  default    : break;
  }
  return false;
}

void MI_Handler::SetScaleMin(double scalemin,unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic:
    p_amisic->HardBase()->SetStop(scalemin,i);
    p_amisic->SoftBase()->SetStart(scalemin,i);
    break;
#endif
  default:
    break;
  }
}

void MI_Handler::SetScaleMax(double scalemax,unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic:
    p_amisic->HardBase()->SetStart(scalemax,i);
    break;
#endif
  default:
    break;
  }
}

double MI_Handler::ScaleMin(unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->HardBase()->Stop(i);
#endif
  default    : break;
  }
  return 0.;
}

double MI_Handler::ScaleMax(unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->HardBase()->Start(i);
#endif
  default    : break;
  }
  return 0.;
}

void MI_Handler::Reset()
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: p_amisic->Reset();
#endif
  default    : break;
  }
}

void MI_Handler::CleanUp()
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: p_amisic->CleanUp();
#endif
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
#ifdef USING__Amisic
  case Amisic: return std::string("Amisic");
#endif
  case None  : return std::string("None");
  default    : break;
  }
  return std::string("Unknown");
}

unsigned int MI_Handler::NIn()
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->HardXS()->NIn();
#endif
  default    : break;
  }
  return 0;
}

unsigned int MI_Handler::NOut()
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->HardXS()->NOut();
#endif
  default    : break;
  }
  return 0;
}

PDF::ISR_Handler *MI_Handler::ISRHandler()
{
  return p_isr;
}

