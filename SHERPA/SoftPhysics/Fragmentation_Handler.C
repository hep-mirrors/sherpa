#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Smart_Pointer.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Exception.H"
#include "AHADIC++/Main/Ahadic.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;
using namespace AHADIC;
#include "AHADIC++/Tools/Hadron_Init.H"

Fragmentation_Handler::Fragmentation_Handler(string _dir,string _file):
  m_dir(_dir), m_file(_file), m_mode(0)
  ,p_ahadic(NULL)
#ifdef USING__PYTHIA
  ,p_lund(NULL)
#endif
{
  Default_Reader reader;
  reader.AddIgnore("[");
  reader.AddIgnore("]");
  reader.SetInputPath(m_dir);
  reader.SetInputFile(m_file);
  m_fragmentationmodel = reader.GetStringNormalisingNoneLikeValues("FRAGMENTATION", string("Ahadic"));
  m_shrink=reader.Get<int>("COMPRESS_PARTONIC_DECAYS",1);
  m_flagpartonics=reader.Get<int>("FLAG_PARTONIC_DECAYS",1);
  if (m_fragmentationmodel==string("Lund")) {
#ifndef USING__PYTHIA
    THROW(fatal_error, "Fragmentation/decay interface to Pythia has not been "+
          string("enabled during compilation (./configure --enable-pythia)."));
#else
    m_sfile=reader.GetValue<string>("LUND_FILE",string("Lund.dat"));
    Hadron_Init init;
    init.Init();
    init.OverWriteProperties(reader);
    ATOOLS::OutputHadrons(msg->Tracking());
    p_lund = new Lund_Interface(m_dir,m_sfile);
    m_mode=1;
    exh->AddTerminatorObject(this);
    // hack for particle initialization, because we don't want to replicate
    // this method in the obsolete Lund Interface.
    return;
#endif
  }
  else if (m_fragmentationmodel==string("Ahadic")) {
    m_sfile=reader.Get<string>("AHADIC_FILE",m_file);
    Hadron_Init init;
    init.Init();
    init.OverWriteProperties(reader);
    ATOOLS::OutputHadrons(msg->Tracking());
    p_ahadic = new AHADIC::Ahadic(m_dir,m_sfile);
    m_mode=2;
    exh->AddTerminatorObject(this);
    return;
  }
  else if (m_fragmentationmodel == string("None")) return;
  THROW(critical_error,"Fragmentation model not implemented.");
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
#ifdef USING__PYTHIA
  if (p_lund!=NULL)   { delete p_lund;   p_lund   = NULL;   }
#endif
  if (p_ahadic!=NULL) { delete p_ahadic; p_ahadic = NULL;   }
  exh->RemoveTerminatorObject(this);
}

void Fragmentation_Handler::PrepareTerminate() 
{
  std::string path(rpa->gen.Variable("SHERPA_STATUS_PATH"));
  if (path=="") return;
  Copy(m_dir+"/"+m_sfile,path+"/"+m_sfile);
}

Return_Value::code 
Fragmentation_Handler::PerformFragmentation(Blob_List *bloblist,
					    Particle_List *particlelist) 
{
  Return_Value::code success;
  switch (m_mode) {
#ifdef USING__PYTHIA
  case 1  : 
    success = p_lund->Hadronize(bloblist);
    if (m_shrink>0 && success==Return_Value::Success) Shrink(bloblist);
    return success;
#endif
  case 2  : 
    success = p_ahadic->Hadronize(bloblist);
    if (success!=Return_Value::Success &&
	success!=Return_Value::Nothing) {
      msg_Debugging()<<"Potential problem in "<<METHOD<<": "<<success<<".\n";
    }
    if (m_shrink>0 && success==Return_Value::Success) Shrink(bloblist);
    return success;
  default : 
    msg_Error()<<"ERROR in "<<METHOD<<":\n"
	       <<"   Unknown hadronization model in mode = "<<m_mode<<".\n"
	       <<"   Abort the run.\n";
    THROW(critical_error,"Fragmentation model not implemented.");
  }
}

void Fragmentation_Handler::Shrink(Blob_List * bloblist) {
  list<Blob *> deleteblobs;
  Particle_Vector * parts;
  for (Blob_List::reverse_iterator blit=bloblist->rbegin();
       blit!=bloblist->rend();++blit) {
    Blob * blob = (*blit);
    if (blob->Type()==btp::Fragmentation) {
      Blob * showerblob(blob->InParticle(0)->ProductionBlob());
      Blob * decblob(showerblob->InParticle(0)->ProductionBlob());
      if (decblob->Type()!=btp::Hadron_Decay) continue;
      showerblob->DeleteInParticles(0);
      showerblob->DeleteOutParticles(0);
      deleteblobs.push_back(blob);
      deleteblobs.push_back(showerblob);
      while (!blob->GetOutParticles().empty()) {
	Particle * part = 
	  blob->RemoveOutParticle(blob->GetOutParticles().front());
	decblob->AddToOutParticles(part);
      }
      decblob->SetStatus(blob_status::needs_hadrondecays);
      decblob->AddData("Partonic",new Blob_Data<int>(m_flagpartonics));
    }
  }
  for (list<Blob *>::iterator blit=deleteblobs.begin();
       blit!=deleteblobs.end();blit++) bloblist->Delete((*blit));
}


