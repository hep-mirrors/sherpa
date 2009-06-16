#include "AMISIC++/Main/Amisic.H"

#include "AMISIC++/Model/Simple_Chain.H"
#include "AMISIC++/Model/Simple_String.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace AMISIC;
using namespace ATOOLS;

Amisic::Amisic():
  m_hardmodel("Unknown"),
  m_softmodel("Unknown"),
  p_hardbase(NULL),
  p_softbase(NULL),
  p_model(NULL),
  p_beam(NULL),
  p_isr(NULL),
  m_external(false) {}

Amisic::Amisic(MODEL::Model_Base *const model,
	       BEAM::Beam_Spectra_Handler *const beam,
	       PDF::ISR_Handler *const isr):
  m_hardmodel("Unknown"),
  m_softmodel("Unknown"),
  p_hardbase(NULL),
  p_softbase(NULL),
  p_model(model),
  p_beam(beam),
  p_isr(isr),
  m_external(true) {}

Amisic::~Amisic() 
{
  if (p_hardbase!=NULL) delete p_hardbase;
  if (p_softbase!=NULL) delete p_softbase;
}

bool Amisic::Initialize()
{
  if (InputPath()=="" && InputFile()=="") return false;
  Data_Reader *reader = new Data_Reader(" ",";","!","=");
  reader->AddWordSeparator("\t");
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile());
  std::vector<std::string> model;
  if (!reader->VectorFromFile(model,"HARD_MODEL_NAME")) {
    model.push_back("Simple_Chain");
  }
  for (size_t i=1;i<model.size();++i) model[0]+=" "+model[i];
  SelectHardModel(model[0]);
  if (!reader->VectorFromFile(model,"SOFT_MODEL_NAME")) {
    model.push_back("None");
  }
  for (size_t i=1;i<model.size();++i) model[0]+=" "+model[i];
  SelectSoftModel(model[0]);
  std::string file;
  if (!reader->ReadFromFile(file,"HARD_MODEL_FILE")) file=InputFile();
  p_hardbase->SetInputPath(InputPath());
  p_hardbase->SetInputFile(file);
  if (!reader->ReadFromFile(file,"SOFT_MODEL_FILE")) file=InputFile();
  p_softbase->SetInputPath(InputPath());
  p_softbase->SetInputFile(file);
  delete reader;
  bool success=true;
  success=success&&p_hardbase->Initialize();
  success=success&&p_softbase->Initialize();
  return success;
}

void Amisic::SameHardProcess(ATOOLS::Blob *blob)
{
  p_hardbase->FillBlob(blob);
}

void Amisic::SameSoftProcess(ATOOLS::Blob *blob)
{
  p_softbase->FillBlob(blob);
}

bool Amisic::VetoHardProcess(ATOOLS::Blob *blob)
{
  return p_hardbase->VetoProcess(blob);
}

bool Amisic::GenerateHardProcess(ATOOLS::Blob *blob)
{
  if (MI_Base::StopGeneration(MI_Base::HardEvent)) return false;
  if (!p_hardbase->DiceProcess()) return false;
  p_hardbase->UpdateAll(p_hardbase);
  return p_hardbase->FillBlob(blob);
}

bool Amisic::GenerateSoftProcess(ATOOLS::Blob *blob)
{
  if (MI_Base::StopGeneration(MI_Base::SoftEvent)) return false;
  if (!p_softbase->DiceProcess()) return false;
  p_softbase->UpdateAll(p_softbase);
  return p_softbase->FillBlob(blob);
}

bool Amisic::GenerateHardEvent(ATOOLS::Blob_List *blobs)
{
  p_hardbase->Reset();
  while (true) {
    Blob *newblob = new Blob();
    if (GenerateHardProcess(newblob)) {
      newblob->SetType(btp::Hard_Collision);
      newblob->SetStatus(blob_status::needs_showers &
			 blob_status::needs_beams &
			 blob_status::needs_hadronization);
      newblob->SetId();
      blobs->push_back(newblob);
    }
    else {
      delete newblob;
      if (MI_Base::StopGeneration(MI_Base::HardEvent)) return true;
      msg_Tracking()<<"Amisic::GenerateHardEvent(): "
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
    Blob *newblob = new Blob();
    if (GenerateSoftProcess(newblob)) {
      newblob->SetType(btp::Soft_Collision);
      newblob->SetStatus(blob_status::needs_beams &
			 blob_status::needs_hadronization);
      newblob->SetId();
      blobs->push_back(newblob);
    }
    else {
      delete newblob;
      if (MI_Base::StopGeneration(MI_Base::SoftEvent)) return true;
      msg_Tracking()<<"Amisic::GenerateSoftEvent(): "
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

void Amisic::CleanUp()
{
  MI_Base::CleanUp();
}

bool Amisic::SelectHardModel(const std::string &hardmodel)
{ 
  m_hardmodel=hardmodel; 
  if (p_hardbase!=NULL) delete p_hardbase;
  msg_Tracking()<<"Amisic::SelectHardModel("<<hardmodel<<"): ";
  if (m_hardmodel=="Simple_Chain") {
    msg_Tracking()<<"Initialize simple hard underlying event model."
		  <<std::endl;
    if (m_external) p_hardbase = new Simple_Chain(p_model,p_beam,p_isr);
    else p_hardbase = new Simple_Chain();
  }
  else {
    msg_Tracking()<<"Initialize no hard underlying event handler."<<std::endl;
    p_hardbase = new MI_None(MI_Base::HardEvent);
    m_hardmodel="None";
  }
  p_hardbase->SetInputPath(InputPath());
  p_hardbase->SetOutputPath(OutputPath());
  p_hardbase->SetInputFile(InputFile());
  return true;
}

bool Amisic::SelectSoftModel(const std::string &softmodel)
{ 
  m_softmodel=softmodel; 
  if (p_softbase!=NULL) delete p_softbase;
  msg_Tracking()<<"Amisic::SelectSoftModel("<<softmodel<<"): ";
  if (m_softmodel=="Simple_String") {
    msg_Tracking()<<"Initialize simple soft underlying event model."
		  <<std::endl;
    if (m_external) p_softbase = new Simple_String();
    else p_softbase = new Simple_String(p_isr);
  }
  else {
    msg_Tracking()<<"Initialize no soft underlying event handler."<<std::endl;
    p_softbase = new MI_None(MI_Base::SoftEvent);
    m_softmodel="None";
  }
  p_softbase->SetInputPath(InputPath());
  p_softbase->SetOutputPath(OutputPath());
  p_softbase->SetInputFile(InputFile());
  return true;
}
