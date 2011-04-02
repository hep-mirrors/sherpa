#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#ifdef USING__PYTHIA
#include "SHERPA/LundTools/Lund_Interface.H"
#endif
#include "ATOOLS/Phys/Mass_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "HADRONS++/Main/Hadrons.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Mixing_Handler.H"


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
  ,p_ampl(NULL)
{
  p_cans = new set<kf_code>;
  Decay_Map* decmap = p_hadrons->DecayMap();
  for (Decay_Map::iterator decit=decmap->begin(); decit!=decmap->end(); decit++)
    p_cans->insert(decit->first.Kfcode());
}

#ifdef USING__PYTHIA
Hadron_Decay_Handler::Hadron_Decay_Handler(Lund_Interface * _lund) :
  m_decmodel(string("Lund")), m_mode(0),
  p_hadrons(NULL), 
  p_lund(_lund), p_ampl(NULL)
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
  // everything prepared that the GenerateMass method with its
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

bool Hadron_Decay_Handler::GenerateMass(ATOOLS::Particle* part, double min, double max) 
{
  DEBUG_FUNC(part->RefFlav()<<" "<<min<<" "<<max);
  double mass = 0.0;
  switch (m_mode) {
#ifdef USING__PYTHIA
  case 0:
    mass = p_lund->GenerateMass(part->RefFlav().Kfcode(),min,max);
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
      mass=hdc->GenerateMass(min, max);
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
      Hadron_Decay_Map* decaymap = p_hadrons->DecayMap();
      Decay_Table* dt = decaymap->FindDecay(blob->InParticle(0)->Flav());
      if(dt->GetDecayChannel(decayproducts)) return true;
      else                                   return false;
    }
  }
  return false;
}

Cluster_Amplitude *Hadron_Decay_Handler::ClusterConfiguration(Blob *const bl)
{
  msg_Debugging()<<METHOD<<"() {\n";
  msg_Indent();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(this);
  for (int i(0);i<bl->NInP();++i) {
    Particle *p(bl->InParticle(i));
    ColorID col(p->GetFlow(2),p->GetFlow(1));
    p_ampl->CreateLeg(-p->Momentum(),p->Flav().Bar(),col,1<<i);
  }
  p_ampl->SetNIn(bl->NInP());
  for (int i(0);i<bl->NOutP();++i) {
    Particle *p(bl->OutParticle(i));
    ColorID col(p->GetFlow(1),p->GetFlow(2));
    p_ampl->CreateLeg(p->Momentum(),p->Flav(),col,1<<(i+p_ampl->NIn()));
  }
  while (p_ampl->Legs().size()>p_ampl->NIn()+2) {
    msg_Debugging()<<*p_ampl<<"\n";
    Cluster_Amplitude *ampl(p_ampl);
    p_ampl = p_ampl->InitNext();
    p_ampl->SetMS(this);
    for (size_t i(0);i<ampl->NIn();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }
    p_ampl->SetNIn(ampl->NIn());
    Cluster_Leg *lij(NULL);
    for (size_t i(ampl->NIn());i<ampl->Legs().size()-1;++i) {
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(i+1);j<ampl->Legs().size();++j) {
	Cluster_Leg *lj(ampl->Leg(j));
	ColorID nc;
	if (li->Col().m_i==0 && li->Col().m_j==0) {
	  nc=lj->Col();
	}
	else if (lj->Col().m_i==0 && lj->Col().m_j==0) {
	  nc=li->Col();
	}
	else if (li->Col().m_i && li->Col().m_i==lj->Col().m_j) {
	  nc.m_i=lj->Col().m_i;
	  nc.m_j=li->Col().m_j;
	}
	else if (li->Col().m_j && li->Col().m_j==lj->Col().m_i) {
	  nc.m_i=li->Col().m_i;
	  nc.m_j=lj->Col().m_j;
	}
	if (nc.m_i>=0 && nc.m_j>=0) {
	  Flavour fl(kf_photon);
	  if (nc.m_i && nc.m_j) fl=Flavour(kf_gluon);
	  else if (nc.m_i) fl=Flavour(kf_d);
	  else if (nc.m_j) fl=Flavour(kf_d).Bar();
	  p_ampl->CreateLeg(li->Mom()+lj->Mom(),fl,nc,li->Id()+lj->Id());
	  lij=p_ampl->Legs().back();
	  break;
	}
      }
      if (lij) break;
    }
    if (lij==NULL) THROW(fatal_error,"Internal eror");
    for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      if (cl->Id()&lij->Id()) continue;
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }    
  }
  double mu2=p_ampl->Leg(0)->Mom().Abs2();
  p_ampl->SetMuF2(mu2);
  p_ampl->SetKT2(mu2);
  msg_Debugging()<<*p_ampl<<"\n";
  while (p_ampl->Prev()) {
    p_ampl=p_ampl->Prev();
    p_ampl->SetMuF2(mu2);
    p_ampl->SetKT2(mu2);
  }
  msg_Debugging()<<"}\n";
  return p_ampl;
}

double Hadron_Decay_Handler::Mass(const ATOOLS::Flavour &fl) const
{
  return fl.HadMass();
}
