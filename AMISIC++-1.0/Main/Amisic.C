#include "Amisic.H"

using namespace AMISIC;

Amisic::Amisic():
  m_hardmodel(Unknown),
  m_softmodel(Unknown),
  p_hardbase(NULL),
  p_softbase(NULL),
  p_model(NULL),
  p_beam(NULL),
  p_isr(NULL),
  m_external(false) {}

Amisic::Amisic(MODEL::Model_Base *_p_model,
	       BEAM::Beam_Spectra_Handler *_p_beam,PDF::ISR_Handler *_p_isr):
  m_hardmodel(Unknown),
  m_softmodel(Unknown),
  p_hardbase(NULL),
  p_softbase(NULL),
  p_model(_p_model),
  p_beam(_p_beam),
  p_isr(_p_isr),
  m_external(true) {}

Amisic::~Amisic() 
{
  if (p_hardbase!=NULL) delete p_hardbase;
  if (p_softbase!=NULL) delete p_softbase;
}

bool Amisic::Initialize()
{
  if (!CheckInputPath()) return false;
  if (!CheckInputFile()) return false;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetFileName(InputFile());
  std::vector<std::string> model;
  if (!reader->VectorFromFile(model,"HARD_MODEL_NAME",ATOOLS::noinputtag,reader->VHorizontal)) {
    model.push_back("None");
  }
  for (unsigned int i=1;i<model.size();++i) model[0]+=std::string(" ")+model[i];
  SelectHardModel(StringToModelID(model[0]));
  if (!reader->VectorFromFile(model,"SOFT_MODEL_NAME",ATOOLS::noinputtag,reader->VHorizontal)) {
    model.push_back("None");
  }
  for (unsigned int i=1;i<model.size();++i) model[0]+=std::string(" ")+model[i];
  SelectSoftModel(StringToModelID(model[0]));
  std::string file;
  reader->ReadFromFile(file,"HARD_MODEL_FILE");
  p_hardbase->SetInputPath(InputPath());
  p_hardbase->SetInputFile(file);
  reader->ReadFromFile(file,"SOFT_MODEL_FILE");  
  p_softbase->SetInputPath(InputPath());
  p_softbase->SetInputFile(file);
  delete reader;
  p_hardbase->Initialize();
  p_softbase->Initialize();
  return true;
}

void Amisic::SameHardProcess(ATOOLS::Blob *blob)
{
  p_hardbase->CreateBlob(blob);
}

void Amisic::SameSoftProcess(ATOOLS::Blob *blob)
{
  p_softbase->CreateBlob(blob);
}

bool Amisic::GenerateHardProcess(ATOOLS::Blob *blob)
{
  if (!p_hardbase->DiceOrderingParameter()) return false;
  if (!p_hardbase->DiceProcess()) return false;
  return p_hardbase->CreateBlob(blob);
}

bool Amisic::GenerateSoftProcess(ATOOLS::Blob *blob)
{
  if (!p_softbase->DiceOrderingParameter()) return false;
  if (!p_softbase->DiceProcess()) return false;
  return p_softbase->CreateBlob(blob);
}

bool Amisic::GenerateHardEvent(ATOOLS::Blob_List *blobs)
{
  p_hardbase->Reset();
  while (true) {
    ATOOLS::Blob *newblob = new ATOOLS::Blob();
    if (GenerateHardProcess(newblob)) {
      newblob->SetId(blobs->size());
      blobs->push_back(newblob);
    }
    else {
      delete newblob;
      if (MI_Base::StopGeneration(MI_Base::HardEvent)) return true;
      ATOOLS::msg.Tracking()<<"Amisic::GenerateHardEvent(): "
			    <<"Cannot create hard underlying event."<<std::endl
			    <<"   Abort attempt."<<std::endl;
      return false;
    }
  }
  return true;
}

bool Amisic::GenerateSoftEvent(ATOOLS::Blob_List *blobs)
{
  p_softbase->Reset();
  while (true) {
    ATOOLS::Blob *newblob = new ATOOLS::Blob();
    if (GenerateSoftProcess(newblob)) {
      newblob->SetId(blobs->size());
      blobs->push_back(newblob);
    }
    else {
      delete newblob;
      if (MI_Base::StopGeneration(MI_Base::SoftEvent)) return true;
      ATOOLS::msg.Tracking()<<"Amisic::GenerateSoftEvent(): "
			    <<"Cannot create soft underlying event."<<std::endl
			    <<"   Abort attempt."<<std::endl;
      return false;
    }
  } 
  return true;
}

bool Amisic::GenerateEvent(ATOOLS::Blob_List *blobs)
{
  Reset();
  if (!GenerateHardEvent(blobs)) return false;
  if (!GenerateSoftEvent(blobs)) return false;
  return true;
}

void Amisic::Reset()
{
  MI_Base::ResetAll();
}

bool Amisic::SelectHardModel(ModelID _m_hardmodel)
{ 
  m_hardmodel=_m_hardmodel; 
  if (p_hardbase!=NULL) delete p_hardbase;
  switch (m_hardmodel) {
  case SimpleChain:
    ATOOLS::msg.Tracking()<<"Amisic::SelectHardModel("<<_m_hardmodel<<"): "
			  <<"Initialize simple hard underlying event model."<<std::endl;
    if (m_external) p_hardbase = new Simple_Chain(p_model,p_beam,p_isr);
    else p_hardbase = new Simple_Chain();
    break;
  case None:
    ATOOLS::msg.Tracking()<<"Amisic::SelectHardModel("<<_m_hardmodel<<"): "
			  <<"Initialize no hard underlying event handler."<<std::endl;
    p_hardbase = new MI_None(MI_Base::HardEvent);
    break;
  default:
    ATOOLS::msg.Error()<<"Amisic::SelectHardModel("<<m_hardmodel<<"): "
		       <<"Model type inconsistent!"<<std::endl
		       <<"   Initialize no hard underlying event handler."<<std::endl;
    m_hardmodel=None;
    p_hardbase = new MI_None(MI_Base::HardEvent);
    break;
  }
  p_hardbase->SetInputPath(m_inputpath);
  p_hardbase->SetOutputPath(m_outputpath);
  p_hardbase->SetInputFile(m_inputfile);
  return true;
}

bool Amisic::SelectSoftModel(ModelID _m_softmodel)
{ 
  m_softmodel=_m_softmodel; 
  if (p_softbase!=NULL) delete p_softbase;
  switch (m_softmodel) {
  case SimpleChain:
    ATOOLS::msg.Tracking()<<"Amisic::SelectSoftModel("<<_m_softmodel<<"): "
			  <<"Initialize simple hard underlying event model."<<std::endl;
    if (m_external) p_softbase = new Simple_Chain(p_model,p_beam,p_isr);
    else p_softbase = new Simple_Chain();
    break;
  case None:
    ATOOLS::msg.Tracking()<<"Amisic::SelectSoftModel("<<_m_softmodel<<"): "
			  <<"Initialize no soft underlying event handler."<<std::endl;
    p_softbase = new MI_None(MI_Base::SoftEvent);
    break;
  default:
    ATOOLS::msg.Error()<<"Amisic::SelectSoftModel("<<m_softmodel<<"): "
		       <<"Model type inconsistent!"<<std::endl
		       <<"   Initialize no soft underlying event handler."<<std::endl;
    m_softmodel=None;
    p_softbase = new MI_None(MI_Base::SoftEvent);
    break;
  }
  p_softbase->SetInputPath(m_inputpath);
  p_softbase->SetOutputPath(m_outputpath);
  p_softbase->SetInputFile(m_inputfile);
  return true;
}

Amisic::ModelID Amisic::StringToModelID(std::string model) 
{
  if (model==std::string("Simple Chain")) return SimpleChain;
  if (model==std::string("None"))         return None;
  ATOOLS::msg.Error()<<"Amisic::StringToModelID("<<model<<"): "
		     <<"Model type unknown!"<<std::endl
		     <<"   Returning 'Unknown'."<<std::endl;
  return Unknown;
}

std::string Amisic::ModelIDToString(ModelID model) 
{
  switch (model) {
  case SimpleChain: return std::string("Simple Chain");
  case None       : return std::string("None");
  case Unknown    : return std::string("Unknown");
  }
  return std::string("Unknown");
}
