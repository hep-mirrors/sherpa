#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"
#include "AHADIC++/Main/Ahadic.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;
using namespace AHADIC;
#include "AHADIC++/Tools/Hadron_Init.H"

Fragmentation_Handler::Fragmentation_Handler(string _shower):
  m_mode(0)
  ,p_ahadic(NULL)
#ifdef USING__PYTHIA
  ,p_lund(NULL)
#endif
{
  Settings& s = Settings::GetMainSettings();
  m_fragmentationmodel = s["FRAGMENTATION"].Get<std::string>();
  m_shrink = s["COMPRESS_PARTONIC_DECAYS"].SetDefault(true).Get<bool>();
  m_flagpartonics = s["FLAG_PARTONIC_DECAYS"].SetDefault(true).Get<bool>();
  if (m_fragmentationmodel==string("Lund")) {
#ifndef USING__PYTHIA
    THROW(fatal_error, "Fragmentation/decay interface to Pythia has not been "+
          string("enabled during compilation (./configure --enable-pythia)."));
#else
    Hadron_Init().Init();
    ATOOLS::OutputHadrons(msg->Tracking());
    p_lund = new Lund_Interface();
    m_mode=1;
    exh->AddTerminatorObject(this);
    // hack for particle initialization, because we don't want to replicate
    // this method in the obsolete Lund Interface.
    return;
#endif
  }
  else if (m_fragmentationmodel==string("Ahadic")) {
    Hadron_Init().Init();
    ATOOLS::OutputHadrons(msg->Tracking());
    p_ahadic = new AHADIC::Ahadic(_shower);
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

Return_Value::code Fragmentation_Handler::
operator()(Blob_List *bloblist, Particle_List *particlelist) 
{   
  Return_Value::code success = Return_Value::Nothing;
  if (m_mode==0 || bloblist->size()==0) return success;
  switch (m_mode) {
#ifdef USING__PYTHIA
  case 1  : 
    success = p_lund->Hadronize(bloblist);
    break;
#endif
  case 2  : 
    success = p_ahadic->Hadronize(bloblist);
    break;
  default : 
    msg_Error()<<"ERROR in "<<METHOD<<":\n"
	       <<"   Unknown hadronization model in mode = "<<m_mode<<".\n"
	       <<"   Abort the run.\n";
    THROW(critical_error,"Fragmentation model not implemented.");
  }
  if (success!=Return_Value::Success &&
      success!=Return_Value::Nothing) {
    msg_Debugging()<<"Potential problem in "<<METHOD<<": "<<success<<".\n";
  }
  if (m_shrink && success==Return_Value::Success) Shrink(bloblist);
  return success;
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
