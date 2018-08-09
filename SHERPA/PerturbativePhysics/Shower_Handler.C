#include "SHERPA/PerturbativePhysics/Shower_Handler.H"

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace PDF;
using namespace ATOOLS;

Shower_Handler::Shower_Handler(const std::string &dir,const std::string &file,
			       MODEL::Model_Base *const model,
			       ISR_Handler *const isr,const int type):
  p_shower(NULL), p_isr(isr)
{
  Default_Reader reader;
  reader.SetInputPath(dir);
  reader.SetInputFile(file);
  m_name=reader.GetStringNormalisingNoneLikeValues("SHOWER_GENERATOR","Dire");
  rpa->gen.SetVariable("JET_CRITERION",
		       reader.Get<std::string>("JET_CRITERION",m_name));
  p_shower =
    Shower_Getter::GetObject(m_name,Shower_Key(model,p_isr,&reader,type));
  if (p_shower==NULL && m_name!="None" &&
      s_loader->LoadLibrary("Sherpa"+m_name)) {
    p_shower =
      Shower_Getter::GetObject(m_name,Shower_Key(model,p_isr,&reader,type));
  }
  if (p_shower==NULL) msg_Info()<<METHOD<<"(): No shower selected."<<std::endl;
}


Shower_Handler::~Shower_Handler() 
{
  if (p_shower) delete p_shower;
}


void Shower_Handler::FillBlobs(ATOOLS::Blob_List * _bloblist) 
{
  if (p_shower && p_shower->ExtractPartons(_bloblist)) {
    Blob * showerblob = _bloblist->FindLast(btp::Shower);
    if (!showerblob->MomentumConserved() && msg->LevelIsDebugging()) {
      msg_Out()<<"Error in "<<METHOD<<": shower violates four-momentum conservation "
	       <<showerblob->CheckMomentumConservation()<<":\n"
	       <<(*showerblob);
    }
    return;
  }
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

void Shower_Handler::SetRemnants(REMNANTS::Remnant_Handler* remnants)
{
  if (p_shower)
    p_shower->SetRemnants(remnants);
}
