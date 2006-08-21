#include "Hadron_Decay_Handler.H"
#include "Message.H"
#include "Random.H"
#include "Blob.H"
#include "Lund_Interface.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef USING__Hadrons
#include "Hadrons.H"
using namespace HADRONS;
#endif


#ifdef USING__Hadrons
Hadron_Decay_Handler::Hadron_Decay_Handler(Hadrons * _hadrons) :
  m_decmodel(string("Hadrons")), m_mode(1),
  p_hadrons(_hadrons),
  p_lund(NULL)
{
  p_cans = new set<kf::code>;
  map<kf::code,Decay_Table *> * decmap = p_hadrons->GetDecayMap();
  for (map<kf::code,Decay_Table *>::iterator decit=decmap->begin();
       decit!=decmap->end();decit++) p_cans->insert(decit->first);
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
  p_lund->SetAllStable();
  p_lund->SwitchOffMassSmearing();
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
    if (p_cans->find(kf)!=p_cans->end()) return true;
    return false;
#ifdef USING__Hadrons
  case 1:
    if (p_cans->find(kf)!=p_cans->end()) return true;
    return false;
#endif
  }
  return false;
}

Return_Value::code Hadron_Decay_Handler::FillHadronDecayBlob(Blob *blob)
{
  Return_Value::code ret = Return_Value::Success;
  switch (m_mode) {
#ifdef USING__Hadrons
    case 1:
      blob->SetTypeSpec("Sherpa");
      ret = p_hadrons->PerformSingleDecay(blob);
      break;
#endif
    case 0:
      blob->SetTypeSpec("Pythia_v6.214");
      ret = p_lund->PerformDecay(blob);
      break;
  }
  return ret;
}
