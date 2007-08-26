#include "Hadron_Decay_Handler.H"
#include "Message.H"
#include "Random.H"
#include "Blob.H"
#include "Particle.H"
#include "Lund_Interface.H"
#include "Mass_Handler.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef USING__Hadrons
#include "Hadrons.H"
#include "Hadron_Decay_Channel.H"
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
      if (p_lund->IsAllowedDecay(flav.Kfcode())) {
        p_cans->insert(flav.Kfcode());
        p_lund->AdjustProperties(flav);
      }
    }
    if( flav.Kfcode()==kf::K_L || flav.Kfcode()==kf::K_S || flav.Kfcode()==kf::K) {
      // adjust for K0, KL and KS even if stable,
      // otherwise 1->1 decay with different masses fails
      p_lund->AdjustProperties(flav);
    }
  }
  p_lund->SwitchOffMassSmearing();
}

Hadron_Decay_Handler::~Hadron_Decay_Handler() 
{
  delete p_cans;
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

Return_Value::code Hadron_Decay_Handler::FillHadronDecayBlob(Blob *blob,const Vec4D& labmom)
{
  Return_Value::code ret = Return_Value::Success;
  switch (m_mode) {
#ifdef USING__Hadrons
    case 1:
      blob->SetTypeSpec("Sherpa");
      ret = p_hadrons->PerformDecay(blob,labmom);
      break;
#endif
    case 0:
      blob->SetTypeSpec("Pythia_v6.214");
      ret = p_lund->PerformDecay(blob);
      break;
  }

  return ret;
}

bool Hadron_Decay_Handler::DiceMass(ATOOLS::Particle* part, double min, double max) {
#ifdef DEBUG__Hadrons
  if(min<0.0 || max<0.0 || !(min<max)) {
    msg_Error()<<METHOD<<" with strange min, max encountered: min="<<min<<" max="<<max<<endl;
  }
#endif
  double mass = 0.0;
  switch (m_mode) {
  case 0:
    mass = p_lund->DiceMass(part->RefFlav().Kfcode(),min,max);
    break;
#ifdef USING__Hadrons
  case 1:
    kf::code kfc = part->RefFlav().Kfcode();
    if(kfc==kf::K || kfc==kf::K_S || kfc==kf::K_L) return true;
    Mass_Handler masshandler(part->RefFlav());
    Blob_Data_Base* data = (*(part->DecayBlob()))["hdc"];
    Hadron_Decay_Channel * hdc=NULL;
    if(data) {
      hdc = data->Get<Hadron_Decay_Channel*>();
      double decaymin = hdc->DecayChannel()->MinimalMass();
      if(decaymin>max) mass = -1.0;
      else             mass = masshandler.GetMass(decaymin, max);
    }
    else mass = masshandler.GetMass(min, max);
    break;
#endif
  }
  
  if(mass>0.0) {
    part->SetFinalMass(mass);
    return true;
  }
  else return false;
}

void Hadron_Decay_Handler::SetSignalProcessBlob(ATOOLS::Blob* spblob)
{
#ifdef USING__Hadrons
  if(m_mode==1) p_hadrons->SetSignalProcessBlob(spblob);
#endif
}
