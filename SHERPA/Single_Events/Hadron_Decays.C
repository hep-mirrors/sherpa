#include "Hadron_Decays.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Mass_Handler.H"
#include "Momenta_Stretcher.H"
#include "Soft_Photon_Handler.H"

#include <utility>
#include <algorithm>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

int Hadron_Decays::m_mass_smearing(1);

Hadron_Decays::Hadron_Decays(HDHandlersMap * _dechandlers, 
                             Soft_Photon_Handler * softphotons) :
  p_dechandlers(_dechandlers), p_softphotons(softphotons), p_bloblist(NULL)
{
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
  msg_Info()<<"Creating mass histograms in masses.root"<<endl;
  p_file = new TFile("masses.root","RECREATE");
  KFCode_ParticleInfo_Map::const_iterator it;
  for (it=s_kftable.begin();it!=s_kftable.end();it++) {
    Flavour flav(it->first);
    if (flav.IsOn() && (flav.IsHadron() || flav.IsLepton())) {
      double limit = flav.Width()==0.0?0.01*flav.PSMass():flav.Width();
      double min = flav.PSMass()-3.0*limit;
      double max = flav.PSMass()+3.0*limit;
      TH1D* myhist = new TH1D(flav.ShellName().c_str(),flav.IDName().c_str(),100,min,max);
      myhist->SetDirectory(p_file);
      mass_hists.insert(make_pair(flav.Kfcode(),myhist));
    }
  }
#endif
#endif
  m_name      = std::string("Hadron_Decays");
  m_type      = eph::Hadronization;
}

Hadron_Decays::~Hadron_Decays()
{
}

Return_Value::code Hadron_Decays::Treat(ATOOLS::Blob_List * bloblist, double & weight) 
{
  DEBUG_FUNC("bloblist->size()="<<bloblist->size());
  if(p_dechandlers->empty() || bloblist->empty()) return Return_Value::Nothing;

  p_bloblist = bloblist;
  for (HDHandlersIter hd=p_dechandlers->begin();hd!=p_dechandlers->end();hd++) {
    hd->second->SetSignalProcessBlob(p_bloblist->FindFirst(btp::Signal_Process));
  }
  bool didit(false);
  for (size_t blit(0);blit<bloblist->size();++blit) {
    Blob* blob=(*bloblist)[blit];
    if (blob->Has(blob_status::needs_hadrondecays)) {
      DEBUG_VAR(*blob);
      if(!blob->MomentumConserved()) {
        msg_Error()<<METHOD<<" found blob with momentum violation: "<<endl
                   <<*blob<<endl
                   <<blob->CheckMomentumConservation()<<endl
                   <<"Will retry event."<<endl;
        return Return_Value::Retry_Event;
      }
      didit = true;
      if(RejectExclusiveChannelsFromFragmentation(blob)) continue;
      
      // random shuffle, against bias in spin correlations and mixing
      Particle_Vector daughters = blob->GetOutParticles();
      random_shuffle(daughters.begin(), daughters.end(), ran);
      
      // fragmentation blobs still contain on-shell particles, stretch them off-shell
      for (Particle_Vector::iterator it=daughters.begin(); it!=daughters.end(); ++it) {
        if((*it)->Status()==part_status::decayed || (*it)->DecayBlob()) {
          msg_Error()<<METHOD<<" Error: Particle "<<*(*it)<<" is already decayed. "
              <<"Will retry event"<<endl;
          return Return_Value::Retry_Event;
        }
        if((*it)->Flav().IsStable()) continue;
        if(!CreateDecayBlob(*it)) return Return_Value::Retry_Event;
      }
      if(!SetMasses(blob)) {
	msg_Error()<<METHOD<<" failed to set masses, retrying event."<<endl<<*blob<<endl;
        return Return_Value::Retry_Event;
      }
      Momenta_Stretcher stretch;
      if(!stretch.StretchBlob(blob)) {
        msg_Error()<<METHOD<<" failed to stretch blob, retrying event."<<endl<<*blob<<endl;
        return Return_Value::Retry_Event;
      }


      for (Particle_Vector::iterator it=daughters.begin(); it!=daughters.end(); ++it) {
        if((*it)->Flav().IsStable()) continue;
        
        if(PerformMixing(*it, p_bloblist)) continue;
        
        if(!Treat((*it)->DecayBlob())) return Return_Value::Retry_Event;
      }
      blob->UnsetStatus(blob_status::needs_hadrondecays);
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
      for (Particle_Vector::iterator it=daughters.begin(); it!=daughters.end(); ++it) {
        if(mass_hists.find((*it)->Flav().Kfcode())!=mass_hists.end()) {
          mass_hists[(*it)->Flav().Kfcode()]->Fill((*it)->Momentum().Mass());
        }
      }
#endif
#endif
    }
  }
  DEBUG_INFO("bloblist->size()="<<bloblist->size()<<" succeeded.");
  return (didit ? Return_Value::Success : Return_Value::Nothing);
}

bool Hadron_Decays::Treat(Blob * blob)
{
  DEBUG_FUNC(*blob);
  Particle* inpart = blob->InParticle(0);
  Vec4D labmom = inpart->Momentum();
  
  // fill decay blob in CMS and all on-shell
  inpart->SetMomentum(Vec4D(inpart->Flav().PSMass(), 0.0, 0.0, 0.0));
  Hadron_Decay_Handler * hdhandler = ChooseDecayHandler(inpart);
  if(hdhandler==NULL || !hdhandler->FillDecayBlob(blob, labmom)) return false;
  inpart->SetStatus(part_status::decayed);
  // random shuffle, against bias in spin correlations and mixing
  Particle_Vector daughters = blob->GetOutParticles();
  random_shuffle(daughters.begin(), daughters.end(), ran);
  
  // create daughter decay blob stubs and then set daughter masses
  if(blob->Has(blob_status::needs_hadrondecays)) {
    for (Particle_Vector::iterator it=daughters.begin();it!=daughters.end();++it) {
      (*it)->SetInfo('D');
      if(!CreateDecayBlob(*it)) return false;
    }
    if(!SetMasses(blob)) return false;
  }
  
  if(!BoostAndStretch(blob, labmom)) return false;
  SetPosition(blob);

  // add extra QED radiation
  if(!AttachExtraQED(blob)) {
    msg_Error()<<METHOD<<": Attaching extra QED failed. Retry event..."
               <<endl;
    return false;
  }

  // recursively treat the created daughter decay blobs
  if(blob->Has(blob_status::needs_hadrondecays)) {
    for (Particle_Vector::iterator it=daughters.begin(); it!=daughters.end(); ++it) {
      if((*it)->Flav().IsStable()) continue;
      if((*it)->Status()==part_status::decayed) {
        msg_Error()<<METHOD<<" Error: Particle "<<*(*it)<<" is already decayed. "
            <<"Will retry event."<<endl;
        return false;
      }
      if(!(*it)->DecayBlob()) {
        msg_Error()<<METHOD<<" Error: Particle "<<*(*it)<<" doesn't have a decay blob stub. "
            <<"Will retry event."<<endl;
        return false;
      }
      if(PerformMixing(*it, p_bloblist)) continue;
      if(!Treat((*it)->DecayBlob())) return false;
    }
  }
  blob->UnsetStatus(blob_status::needs_hadrondecays);
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
  for (Particle_Vector::iterator it=daughters.begin(); it!=daughters.end(); ++it) {
    if(mass_hists.find((*it)->Flav().Kfcode())!=mass_hists.end()) {
      mass_hists[(*it)->Flav().Kfcode()]->Fill((*it)->Momentum().Mass());
    }
  }
#endif
#endif
  DEBUG_INFO("succeeded.");
  return true;
}

bool Hadron_Decays::CreateDecayBlob(Particle* inpart)
{
  DEBUG_FUNC(inpart->Flav());
  if(inpart->DecayBlob()) abort();
  if(inpart->Flav().IsStable()) return true;
  if(inpart->Time()==0.0) inpart->SetTime();
  Blob* blob = p_bloblist->AddBlob(btp::Hadron_Decay);
  blob->AddToInParticles(inpart);
  SetPosition(blob);
  Hadron_Decay_Handler * hdhandler = ChooseDecayHandler(inpart);
  if(hdhandler==NULL || !hdhandler->CreateDecayBlob(blob)) {
    msg_Error()<<"Failed to create decay blob for "<<*inpart<<endl
               <<"Is a decay handler active for this particle?"<<endl<<endl;
    return false;
  }
  blob->SetStatus(blob_status::needs_hadrondecays);
  DEBUG_INFO("succeeded.");
  return true;
}

bool SortByWidth(Particle* p1, Particle* p2) {
  return p1->Flav().Width() < p2->Flav().Width();
}

bool Hadron_Decays::SetMasses(Blob * blob)
{
  DEBUG_FUNC(blob->Id());
  Particle_Vector daughters = blob->GetOutParticles();
  if(m_mass_smearing==0 || daughters.size()==1) return true;
  sort(daughters.begin(), daughters.end(), SortByWidth);
  double max_mass;
  if(blob->NInP()==1) max_mass = blob->InParticle(0)->FinalMass();
  else {
    Vec4D total(0.0,0.0,0.0,0.0);
    for(int i=0; i<blob->NInP(); i++) total+=blob->InParticle(i)->Momentum();
    max_mass = total.Mass();
  }
  
  bool success=true; size_t cnt=0;
  do {
    success=true;
    double max = max_mass;
    for(Particle_Vector::iterator it=daughters.begin();it!=daughters.end();++it) {
      if(m_mass_smearing==2 && (*it)->Flav().IsStable()) continue;
      if(!(*it)->Flav().IsHadron() || (*it)->Flav().IsStable()) {
        Mass_Handler masshandler((*it)->RefFlav());
        double mass = masshandler.GetMass(0.0, max);
        max-=mass;
        (*it)->SetFinalMass(mass);
      }
      else {
        Hadron_Decay_Handler * hdhandler = ChooseDecayHandler((*it));
        if(!hdhandler->DiceMass(*it,0.,max)) {
          success=false; ++cnt;
	  if(cnt>22) { DEBUG_INFO("failed."); return false;}
          break;
        }
        max-=(*it)->FinalMass();
      }
    }
  } while(success==false);

  DEBUG_INFO("succeeded.");
  return true;
}

bool Hadron_Decays::BoostAndStretch(Blob* blob, const Vec4D& labmom)
{
  DEBUG_FUNC(blob->Id()<<", "<<labmom);
  // 1.
  Particle* inpart = blob->InParticle(0);
  double factor = inpart->FinalMass()/inpart->Flav().PSMass();
  inpart->SetMomentum(Vec4D(inpart->FinalMass(), 0., 0., 0.));
  Particle_Vector daughters = blob->GetOutParticles();
  for(Particle_Vector::iterator it=daughters.begin();it!=daughters.end();++it) {
    Vec4D mom = (*it)->Momentum();
    (*it)->SetMomentum(Vec4D(factor*mom[0], Vec3D(mom)));
  }
  
  // 2.
  Momenta_Stretcher stretch;
  if(!stretch.StretchBlob(blob)) {
    msg_Error()<<METHOD<<" failed to stretch blob, retrying event."<<endl<<*blob<<endl;
    return false;
  }

  // 3.
  Poincare boost(labmom);
  boost.Invert();
  blob->Boost(boost);
  
  DEBUG_INFO("succeeded.");
  return true;
}

Hadron_Decay_Handler* Hadron_Decays::ChooseDecayHandler(Particle* part)
{
  for (HDHandlersIter hd=p_dechandlers->begin();hd!=p_dechandlers->end();hd++) {
    if (hd->second->CanDealWith(part->Flav().Kfcode())) {
      return hd->second;
    }
  }
  return NULL;
}

bool Hadron_Decays::RejectExclusiveChannelsFromFragmentation(Blob* fragmentationblob)
{
  if(fragmentationblob->Type()==btp::Fragmentation) {
    Blob* showerblob = fragmentationblob->UpstreamBlob();
    if(showerblob && showerblob->Type()==btp::FS_Shower) {
      Blob* decayblob = showerblob->UpstreamBlob();
      if(decayblob && decayblob->Type()==btp::Hadron_Decay) {
        DEBUG_FUNC(fragmentationblob->Id());
        rvalue.IncCall(METHOD);
        FlavourSet decayresults;
        for(int i=0;i<fragmentationblob->NOutP();i++) {
          decayresults.insert(fragmentationblob->OutParticle(i)->Flav());
        }
        Hadron_Decay_Handler * hdhandler =
            ChooseDecayHandler(decayblob->InParticle(0));
        if(hdhandler->IsExclusiveDecaychannel(decayblob,decayresults)) {
          DEBUG_INFO("found exclusive decay channel, retrying.");
          rvalue.IncRetryPhase(METHOD);
          p_bloblist->Delete(fragmentationblob);
          p_bloblist->Delete(showerblob);
          decayblob->SetStatus(blob_status::needs_hadrondecays);
          decayblob->AddStatus(blob_status::internal_flag);
          decayblob->DeleteOutParticles();
          decayblob->InParticle(0)->SetStatus(part_status::active);
          Treat(decayblob);

          return true;
        }
        DEBUG_INFO("did not find exclusive decay channel, continue.");
      }
    }
  }
  return false;
}

bool Hadron_Decays::PerformMixing(Particle* inpart, Blob_List* bloblist)
{
  Hadron_Decay_Handler * hdhandler = ChooseDecayHandler(inpart);
  if(hdhandler==NULL) return false;
  return hdhandler->PerformMixing(inpart, bloblist);
}

bool Hadron_Decays::AttachExtraQED(Blob* blob)
{
  DEBUG_FUNC(blob->Id());
  // attach QED radiation to blobs before they are subsequently decayed
  if (blob->Status()!=blob_status::needs_hadrondecays) return true;
  if (!p_softphotons) return false;
  blob->SetStatus(blob_status::needs_extraQED);
  bool ret = p_softphotons->AddRadiation(blob);
  blob->SetStatus(blob_status::needs_hadrondecays);
  DEBUG_INFO("succeeded.");
  return ret;
}

void Hadron_Decays::SetPosition(ATOOLS::Blob* blob)
{
  Particle* inpart = blob->InParticle(0);
  if(inpart->Flav().Kfcode()==kf_K) return;
  
  // boost lifetime into lab
  double gamma = 1./rpa.gen.Accu();
  if (inpart->Flav().PSMass()>rpa.gen.Accu()) {
    gamma = inpart->E()/inpart->Flav().PSMass();
  }
  else {
    double q2    = dabs(inpart->Momentum().Abs2());
    if (q2>rpa.gen.Accu()) gamma = inpart->E()/sqrt(q2);
  }
  double lifetime_boosted = gamma * inpart->Time();
  
  Vec3D      spatial = inpart->Distance( lifetime_boosted );
  Vec4D     position = Vec4D( lifetime_boosted*rpa.c(), spatial );
  blob->SetPosition( inpart->XProd() + position ); // in mm
}

void Hadron_Decays::CleanUp()
{
  for (HDHandlersIter hd=p_dechandlers->begin();hd!=p_dechandlers->end();hd++) {
    hd->second->CleanUp();
  }
}

void Hadron_Decays::Finish(const std::string &) 
{
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
  map<kf_code,TH1D*>::iterator it;
  for(it=mass_hists.begin();it!=mass_hists.end();it++) {
    if(it->second->GetEntries()>0) {
      it->second->Write();
    }
  }
  p_file->Close();
#endif
#endif
}
