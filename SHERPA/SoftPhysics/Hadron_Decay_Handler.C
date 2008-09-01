#include "Hadron_Decay_Handler.H"
#include "Message.H"
#include "Random.H"
#include "Blob.H"
#include "Blob_List.H"
#include "Particle.H"
#include "CXXFLAGS.H"
#ifdef USING__PYTHIA
#include "Lund_Interface.H"
#endif
#include "Mass_Handler.H"
#include "Run_Parameter.H"
#include "Hadrons.H"
#include "Decay_Map.H"
#include "Hadron_Decay_Table.H"
#include "Hadron_Decay_Channel.H"
#include "Mixing_Handler.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;
using namespace HADRONS;


Hadron_Decay_Handler::Hadron_Decay_Handler(Hadrons * _hadrons) :
  m_decmodel(string("Hadrons")), m_mode(1),
  p_hadrons(_hadrons)
#ifdef USING__PYTHIA
  ,p_lund(NULL)
#endif
{
  p_cans = new set<kf_code>;
  FlDtVMap* decmap = p_hadrons->DecayMap();
  for (FlDtVMap::iterator decit=decmap->begin(); decit!=decmap->end(); decit++)
    p_cans->insert(decit->first.Kfcode());
}

#ifdef USING__PYTHIA
Hadron_Decay_Handler::Hadron_Decay_Handler(Lund_Interface * _lund) :
  m_decmodel(string("Lund")), m_mode(0),
  p_hadrons(NULL), 
  p_lund(_lund)
{ 
  p_cans = new set<kf_code>;
  Flavour flav(kf_tau);
  if (flav.IsOn() && !flav.IsStable()) {
    if (p_lund->IsAllowedDecay(flav.Kfcode())) p_cans->insert(flav.Kfcode());
  }
  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end();++kfit) {
    Flavour flav(kfit->first);
    if (flav.IsOn() && flav.IsHadron() && !flav.IsStable()) {
      if (p_lund->IsAllowedDecay(flav.Kfcode())) {
        p_cans->insert(flav.Kfcode());
        p_lund->AdjustProperties(flav);
      }
    }
    if( flav.Kfcode()==kf_K_L || flav.Kfcode()==kf_K_S || flav.Kfcode()==kf_K) {
      // adjust for K0, KL and KS even if stable,
      // otherwise 1->1 decay with different masses fails
      p_lund->AdjustProperties(flav);
    }
  }
  p_lund->SwitchOffMassSmearing();
}
#endif

Hadron_Decay_Handler::~Hadron_Decay_Handler() 
{
  delete p_cans;
  if (p_hadrons) delete p_hadrons; p_hadrons=NULL;
}

bool Hadron_Decay_Handler::CanDealWith(kf_code kf) {
  switch (m_mode) {
  case 0:
    if (p_cans->find(kf)!=p_cans->end()) return true;
    return false;
  case 1:
    if (p_cans->find(kf)!=p_cans->end()) return true;
    return false;
  }
  return false;
}

bool Hadron_Decay_Handler::CreateDecayBlob(Blob* blob)
{
  DEBUG_FUNC("blob->Id()="<<blob->Id());
  // after this method has run, the blob is supposed to have 
  // everything prepared that the DiceMass method with its
  // InParticle needs.
  
//   if(part->Time()==0.0) part->SetTime();
//   SetPosition(blob);
  
//   PerformMixing(blob, bloblist);
//   blob=bloblist->back();
  
  bool returncode;
  switch (m_mode) {
    case 1:
      DEBUG_INFO("with Sherpa.");
      blob->SetTypeSpec("Sherpa");
      returncode = p_hadrons->CreateDecayBlob(blob);
      break;
#ifdef USING__PYTHIA
    case 0:
      DEBUG_INFO("with Pythia.");
      blob->SetTypeSpec("Pythia_v6.214");
      returncode = true;
      break;
#endif
  }
  return true;
}

bool Hadron_Decay_Handler::FillDecayBlob(Blob *blob, const Vec4D& labmom)
{
  DEBUG_FUNC("blob->Id()="<<blob->Id());
  // after this method has run, the blob is supposed to be complete
  // with kinematics in CMS, and with on-shell particles.
  switch (m_mode) {
    case 1:
      DEBUG_INFO("with Sherpa.");
      return p_hadrons->FillDecayBlob(blob, labmom);
#ifdef USING__PYTHIA
    case 0:
      DEBUG_INFO("with Pythia.");
      return p_lund->PerformDecay(blob);
#endif
  }
  return false;
}

bool Hadron_Decay_Handler::DiceMass(ATOOLS::Particle* part, double min, double max) 
{
  DEBUG_FUNC(part->RefFlav()<<" "<<min<<" "<<max);
  double mass = 0.0;
  switch (m_mode) {
#ifdef USING__PYTHIA
  case 0:
    mass = p_lund->DiceMass(part->RefFlav().Kfcode(),min,max);
    break;
#endif
  case 1:
    Blob* decayblob=part->DecayBlob();
    kf_code kfc = part->RefFlav().Kfcode();
    if(kfc==kf_K || kfc==kf_K_S || kfc==kf_K_L || decayblob->Type()!=btp::Hadron_Decay) 
      return true;
    Blob_Data_Base* data = (*decayblob)["hdc"];
    if(data) {
      Hadron_Decay_Channel* hdc = data->Get<Hadron_Decay_Channel*>();
      mass=hdc->DiceMass(min, max);
    }
    else {
      Mass_Handler masshandler(part->RefFlav());
      mass = masshandler.GetMass(min, max);
    }
    break;
  }
  
  DEBUG_VAR(mass);
  if(mass>0.0) {
    part->SetFinalMass(mass);
    return true;
  }
  else return false;
}

void Hadron_Decay_Handler::SetSignalProcessBlob(ATOOLS::Blob* spblob)
{
  if(m_mode==1) p_hadrons->SetSignalProcessBlob(spblob);
}

bool Hadron_Decay_Handler::PerformMixing(Particle* inpart, Blob_List* bloblist)
{
  if(m_mode==1) return p_hadrons->MixingHandler()->PerformMixing(inpart, bloblist);
  else          return false;
}

void Hadron_Decay_Handler::CleanUp()
{
  if(m_mode==1) p_hadrons->CleanUp();
}

bool Hadron_Decay_Handler::IsExclusiveDecaychannel(Blob* blob, FlavourSet decayproducts)
{
  if(m_mode==1) {
    if(blob->TypeSpec()=="Sherpa") {
      Decay_Map* decaymap = p_hadrons->DecayMap();
      Hadron_Decay_Table* dt = decaymap->FindDecay(blob->InParticle(0)->Flav());
      if(dt->GetDecayChannel(decayproducts)) return true;
      else                                   return false;
    }
  }
  return false;
}
