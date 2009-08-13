#include "SHERPA/PerturbativePhysics/Shower_Handler.H"

#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Shower_Handler::Shower_Handler
(const std::string &dir,const std::string &file,
 MODEL::Model_Base *const model,PDF::ISR_Handler *const isr):
  p_shower(NULL), p_isr(isr)
{
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(dir);
  dataread.SetInputFile(file);
  m_name=dataread.GetValue<std::string>("SHOWER_GENERATOR","CSS");
  p_shower = PDF::Shower_Getter::GetObject
    (m_name,PDF::Shower_Key(model,p_isr,&dataread));
  if (p_shower==NULL) msg_Info()<<METHOD<<"(): No shower selected."<<std::endl;
}


Shower_Handler::~Shower_Handler() 
{
  if (p_shower) delete p_shower;
}


void Shower_Handler::FillBlobs(ATOOLS::Blob_List * _bloblist) 
{
  if (p_shower && p_shower->ExtractPartons(_bloblist)) return;
  THROW(fatal_error,"Internal error");
}

void Shower_Handler::FillDecayBlobs(ATOOLS::Blob_List * _bloblist) 
{
  if (p_shower && p_shower->ExtractPartons(_bloblist)) return;
  THROW(fatal_error,"Internal error");
}

void Shower_Handler::CleanUp() 
{
  if (p_shower) p_shower->CleanUp();
}

