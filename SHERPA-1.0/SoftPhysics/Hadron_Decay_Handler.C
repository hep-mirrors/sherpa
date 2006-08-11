#include "Hadron_Decay_Handler.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef USING__Hadrons
using namespace HADRONS;
#endif


#ifdef USING__Hadrons
Hadron_Decay_Handler::Hadron_Decay_Handler(Hadrons * _hadrons) :
  m_decmodel(string("Hadrons")), m_mode(1),
  p_hadrons(_hadrons),
  p_lund(NULL)
{
  p_cans = new set<kf::code>;
  //p_hadrons->FillAllowedDecays(p_cans);
  map<kf::code,Decay_Table *> * decmap = p_hadrons->GetDecayMap();
  for (map<kf::code,Decay_Table *>::iterator decit=decmap->begin();
       decit!=decmap->end();decit++) p_cans->insert(decit->first);
  
  SwitchOfLundDecays();
}
#endif


Hadron_Decay_Handler::Hadron_Decay_Handler(Lund_Interface * _lund) :
  m_decmodel(string("Lund")), m_mode(0),
#ifdef USING__Hadrons
  p_hadrons(NULL), 
#endif
  p_lund(_lund)
{ 
  p_cans = new set<kf::code>;
  Flavour flav(kf::tau);
  if (flav.IsOn() && !flav.IsStable()) {
    if (p_lund->IsAllowedDecay(flav.Kfcode())) p_cans->insert(flav.Kfcode());
  }
  Fl_Iter fli;
  for (flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.IsOn() && flav.IsHadron() && !flav.IsStable()) {
      if (p_lund->IsAllowedDecay(flav.Kfcode())) p_cans->insert(flav.Kfcode());
    }
  }
  // To do: Check for interference with hadrons.
  //  for (set<kf::code>::iterator cit=p_cans->begin();cit!=p_cans->end();cit++) {
  //    std::cout<<"Lund can deal with :"<<Flavour((*cit))<<std::endl;
  //}
}

Hadron_Decay_Handler::~Hadron_Decay_Handler() 
{
}

void Hadron_Decay_Handler::EraseTreated(std::set<int> * hadrons)
{
  if (m_mode==0) hadrons->clear();
#ifdef USING__Hadrons
  if (m_mode==1) {
    
    for (set<kf::code>::iterator citer=p_cans->begin();citer!=p_cans->end();citer++) {
      //msg.Debugging()<<"Killing flavours: "
      //<<citer->first<<" ("<<cans->size()<<" ) "<<hadrons->size()<<endl;
      hadrons->erase(int(*citer));
      //msg.Debugging()<<"                  "
      //<<citer->first<<" ("<<cans->size()<<" ) "<<hadrons->size()<<endl;
    }
  }
#endif
}

bool Hadron_Decay_Handler::CanDealWith(kf::code kf) {
  switch (m_mode) {
  case 0:
    return true;
#ifdef USING__Hadrons
  case 1:
    if (p_cans->find(kf)!=p_cans->end()) return true;
    return false;
#endif
  }
  return false;
}


// with spin correlation 
bool Hadron_Decay_Handler::FillHadronDecayBlob(Particle *part,Blob_List *blob_list)
{
  Blob * blob(new Blob());
  blob->SetType(btp::Hadron_Decay);
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetId();
  blob->SetPosition(part->ProductionBlob()->Position());
  blob->AddToInParticles(part);

  switch (m_mode) {
#ifdef USING__Hadrons
    case 1: break;
#endif
    case 0:
      p_lund->PerformDecay(blob);
      break;
  }
  /*
    To do: PerformDecay doesn't work out: Erase blob, keep particle -> Related to smearing.
  */
  blob_list->push_back(blob);
  return true;
}

void Hadron_Decay_Handler::SwitchOfLundDecays()
{
#ifdef USING__Hadrons
  std::map<ATOOLS::kf::code,ATOOLS::Decay_Table *>::iterator dtiter;
  for (dtiter=p_hadrons->GetDecayMap()->begin();
       dtiter!=p_hadrons->GetDecayMap()->end();dtiter++) {
    p_lund->SwitchOfDecays(dtiter->first);
  }
#endif
}
