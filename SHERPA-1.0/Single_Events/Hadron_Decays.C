#include "Hadron_Decays.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Mass_Handler.H"
#include "Momenta_Stretcher.H"
#include "Amplitude_Tensor.H"

#include <utility>
#include <algorithm>

#ifdef PROFILE__Hadron_Decays
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadron_Decays::Hadron_Decays(HDHandlersMap * _dechandlers) :
  p_dechandlers(_dechandlers), p_bloblist(NULL), p_saved_amplitudes(NULL)
{
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
  p_file = new TFile("masses.root","RECREATE");
  Fl_Iter fli;
  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
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
  if(p_dechandlers->size()) m_mass_smearing=p_dechandlers->begin()->second->GetMassSmearing();
  else                      m_mass_smearing=true;
  m_name      = std::string("Hadron_Decays");
  m_type      = eph::Hadronization;
}

Hadron_Decays::~Hadron_Decays()
{
}

Return_Value::code Hadron_Decays::Treat(ATOOLS::Blob_List * bloblist, double & weight) 
{
  msg_Tracking()<<"--------- "<<METHOD<<" - START --------------"<<endl;
  PROFILE_HERE;
  if(p_dechandlers->empty()) return Return_Value::Nothing;

  if (bloblist->empty()) {
    msg_Error()<<"Potential error in "<<METHOD<<endl
      <<"   Incoming blob list contains "<<bloblist->size()<<" entries."<<endl
      <<"   Continue and hope for the best."<<endl;
    return Return_Value::Nothing;
  }

  p_bloblist = bloblist;
  for (HDHandlersIter hd=p_dechandlers->begin();hd!=p_dechandlers->end();hd++) {
    hd->second->SetSignalProcessBlob(p_bloblist->FindFirst(btp::Signal_Process));
  }
  bool found(true), didit(false);
  while (found) {
    found = false;
    int count = 0;
    for (size_t blit(0);blit<bloblist->size();++blit) {
      count++;
      if ((*bloblist)[blit]->Has(blob_status::needs_hadrondecays)) {
        Return_Value::code blob_ret = Treat((*bloblist)[blit]);
        if( blob_ret == Return_Value::Retry_Event ) return Return_Value::Retry_Event;
        didit = true;
        found = true;
      }
    }
  }
  msg_Tracking()<<"--------- Hadron_Decays::Treat - FINISH -------------"<<endl;
  return (didit ? Return_Value::Success : Return_Value::Nothing);
}


Return_Value::code Hadron_Decays::Treat(Blob * blob)
{
#ifdef DEBUG__Hadrons
  if(!blob->MomentumConserved()) {
    PRINT_INFO(" beware: from the beginning momentum not conserved in "
              <<"blob ["<<blob->Id()<<"] of event "<<rpa.gen.NumberOfDicedEvents()<<endl
              <<blob->CheckMomentumConservation()<<endl
              <<(*blob));
  }
  Vec4D vtotal(0.0,0.0,0.0,0.0); double mtotal=0.0;
  for(int i=0;i<blob->NOutP();i++) {
    vtotal += blob->OutParticle(i)->Momentum();
    mtotal += blob->OutParticle(i)->FinalMass();
  }
  if(vtotal.Mass()+1.e-6<mtotal) {
    std::cout.precision(8);
    PRINT_INFO(" beware: mtotal="<<mtotal<<" > vtotal.Mass()="<<vtotal.Mass()<<" in "
             <<"blob ["<<blob->Id()<<"] of event "<<rpa.gen.NumberOfDicedEvents()<<endl
             <<(*blob));
  }
#endif
  m_daughters = blob->GetOutParticles();
  
  if(!PrepareMassSmearing(blob)) return Return_Value::Retry_Event;

  // Decaying all daughters one by one (in CMS*)
  bool retry_all = false;
  int all_again_trials = 0;
  do { // retry mass smearing of all daughters in case a single-mass retry didn't succeed
    retry_all = false;
    for (size_t i=0;i<m_daughters.size();i++) {
      Return_Value::code particle_ret = Treat(i);
      if( particle_ret  == Return_Value::Failure ) {
        ResetMotherAmplitudes(blob);
        all_again_trials++;
        retry_all = true;
        break;
      }
      else if( particle_ret == Return_Value::Retry_Event ) {
        return Return_Value::Retry_Event;
      }
    }
    if( all_again_trials > 10 ) {
      msg_Error()<<"Warning in "<<METHOD<<" in "
        <<"   blob ["<<blob->Id()<<"] of event "<<rpa.gen.NumberOfDicedEvents()<<endl
        <<"   retried mass dicing too often and didn't succeed. "
        <<"   Will retry event."<<endl;
      return Return_Value::Retry_Event;
    }
  } while( retry_all == true ) ;
  // by here all daughters should be on their new mass shells and have either
  // a decayblob or momentum, in CMS.

  // now we stretch the saved_momenta to the FinalMasses of m_daughters
  if(m_mass_smearing && blob->NOutP()>1) {
    Momenta_Stretcher stretch;
    stretch.StretchMomenta( m_daughters, m_saved_momenta);
  }

  // now we boost all decayblobs back to these stretched target momenta
  for (size_t i=0;i<m_daughters.size();i++) {
    if(m_daughters[i]->DecayBlob()) {
      Poincare boost(m_saved_momenta[i]);
      boost.Invert();
      m_daughters[i]->DecayBlob()->Boost(boost);
      m_daughters[i]->SetStatus(part_status::decayed);
#ifdef DEBUG__Hadrons // fill the mass histogram of this flavour
#ifdef USING__ROOT
      if(mass_hists.find(m_daughters[i]->Flav().Kfcode())!=mass_hists.end()) {
        mass_hists[m_daughters[i]->Flav().Kfcode()]->Fill(m_daughters[i]->FinalMass());
      }
#endif
#endif
    }
    else {
      m_daughters[i]->SetMomentum(m_saved_momenta[i]);
#ifdef DEBUG__Hadrons // fill the mass histogram of this flavour
#ifdef USING__ROOT
      if(mass_hists.find(m_daughters[i]->Flav().Kfcode())!=mass_hists.end()) {
        mass_hists[m_daughters[i]->Flav().Kfcode()]->Fill(m_daughters[i]->FinalMass());
      }
#endif
#endif
    }
  }

  // finally set the position and lifetime in case it's a hadron decay blob
  if(blob->Type() == btp::Hadron_Decay) {
    double time        = blob->InParticle(0)->LifeTime(); // in s
    Vec3D      spatial = blob->InParticle(0)->Distance( time );
    Vec4D     position = Vec4D( time*rpa.c(), spatial );
    blob->SetPosition( blob->InParticle(0)->XProd() + position ); // in mm
  }
  blob->UnsetStatus(blob_status::needs_hadrondecays);
  
#ifdef DEBUG__Hadrons
  if(!blob->MomentumConserved()) {
    PRINT_INFO(" beware: from the beginning momentum not conserved in "
              <<"blob ["<<blob->Id()<<"] of event "<<rpa.gen.NumberOfDicedEvents()<<endl
              <<blob->CheckMomentumConservation()<<endl
              <<(*blob));
  }
#endif
  return Return_Value::Success;
}


Return_Value::code Hadron_Decays::Treat(int i)
{
  double max_mass = m_max_mass;
  for( int j=0; j<i; j++) max_mass -= m_daughters[j]->FinalMass();
  if( !(max_mass > 0.0) ) {
#ifdef DEBUG__Hadrons
    PRINT_INFO(" max_mass="<<max_mass<<" for particle "<<m_daughters[i]->Flav()
               <<" < 0. Retry all particles.");
#endif
    return Return_Value::Failure;
  }
  Particle * part = m_daughters[i];
  Mass_Handler masshandler(part->Flav());
  if( part->Status()==part_status::active && !part->Flav().IsStable() ) {
    Hadron_Decay_Handler * hdhandler = ChooseDecayHandler(part);
    if (hdhandler) {
      Return_Value::code ret = Return_Value::Undefined;
      for(int trials=0; trials<100; trials++) {
        Blob* blob;
        // check if particle has a decayblob already, where its
        // decay channel could have been stored in a previous try (HADRONS)
        if(part->DecayBlob()) {
          blob = part->DecayBlob();
          blob->SetStatus(blob_status::needs_hadrondecays);
          part->SetStatus(part_status::active);
          blob->DeleteOutParticles();
        }
        else {
          blob = new Blob();
          blob->SetId();
          blob->SetType(btp::Hadron_Decay);
          blob->SetStatus(blob_status::needs_hadrondecays);
          blob->AddToInParticles(part);
          p_bloblist->push_back(blob);
        }
        // dice FinalMass of part (but only if mass smearing is on)
        if(m_mass_smearing && m_daughters.size()>1) {
          if(! hdhandler->DiceMass(part,0.,max_mass)) {
            for(size_t k=0;k<m_daughters.size();k++) {
              m_daughters[k]->SetStatus(part_status::active);
              if(m_daughters[k]->DecayBlob()) p_bloblist->Delete(m_daughters[k]->DecayBlob());
            }
            return Return_Value::Failure;
          }
        }
        // do the decay in CMS
        part->SetMomentum(Vec4D(part->FinalMass(),0.0,0.0,0.0));
        ret = hdhandler->FillHadronDecayBlob(blob,m_saved_momenta[i]);
        if( ret == Return_Value::Success ) { return Return_Value::Success; }
        else if (ret == Return_Value::Nothing) {
          p_bloblist->Delete(blob);
          part->SetStatus(part_status::active);
          return Return_Value::Success;
        }
        else if ( ret == Return_Value::Error ) {
          msg_Error()<<"Error in "<<METHOD<<":"<<endl
            <<"   Hadron_Decay_Handler "<<hdhandler->Name()<<" failed to decay "<<endl
            <<"   "<<(*part)<<","<<endl
            <<"   it returned "<<ret<<". Will retry event."<<endl;
          return Return_Value::Retry_Event;
        }
      }
      // if after 100 trials still no success: send back to retry all masses
      for(size_t k=0;k<m_daughters.size();k++) {
        m_daughters[k]->SetStatus(part_status::active);
        if(m_daughters[k]->DecayBlob()) p_bloblist->Delete(m_daughters[k]->DecayBlob());
      }
#ifdef DEBUG__Hadrons
      PRINT_INFO("Tried particle "<<m_daughters[i]->Flav()<<" 100 times with max_mass="
                 <<max_mass<<". Retry all particles.");
#endif
      return Return_Value::Failure;
    }
    // if no decay handler found (-> particle quasi stable)
    else
      if(m_mass_smearing && m_daughters.size()>1)
        part->SetFinalMass(masshandler.GetMass(0,max_mass));
  }
  // if particle stable
  else
    if(m_mass_smearing && m_daughters.size()>1)
      part->SetFinalMass(masshandler.GetMass(0,max_mass));
  return Return_Value::Success;
}

bool SortByWidth(Particle* p1, Particle* p2) {
  return p1->Flav().Width() < p2->Flav().Width();
}

bool Hadron_Decays::PrepareMassSmearing(Blob* blob)
{
  // Save amplitude tensor of mother production blob for case of retry
  Blob* motherblob = blob->InParticle(0)->ProductionBlob();
  Blob_Data_Base* scdata = (*motherblob)["amps"];
  if(scdata) {
    if(p_saved_amplitudes) delete p_saved_amplitudes;
    p_saved_amplitudes = new Amplitude_Tensor(*(scdata->Get<Amplitude_Tensor*>()));
  }
  // Sort daughters by width, necessary for mass treatment
  sort(m_daughters.begin(), m_daughters.end(), SortByWidth);
  m_saved_momenta.clear();
  for(size_t i=0;i<m_daughters.size();i++) {
    m_saved_momenta.push_back(m_daughters[i]->Momentum());
  }
  Vec4D total(0.0,0.0,0.0,0.0);
  for(int i=0;i<blob->NOutP();i++) total += blob->OutParticle(i)->Momentum();
  m_max_mass = total.Mass();
  if( !(m_max_mass > 0.0) ) {
    msg_Error()<<METHOD<<" Error: total outgoing mass in this blob:"<<endl
      <<"  m_max_mass="<<m_max_mass<<endl
      <<"  The blob was:"<<endl
      <<(*blob)<<endl
      <<"  and its totalmomentum.Mass2()="<<total.Abs2()<<endl
      <<"  This should not happen, retrying event."<<endl;
    return false;
  }
  return true;
}

void Hadron_Decays::ResetMotherAmplitudes(Blob* blob)
{
  Blob* motherblob = blob->InParticle(0)->ProductionBlob();
  Blob_Data_Base* scdata = (*motherblob)["amps"];
  if(scdata) {
    *(scdata->Get<Amplitude_Tensor*>())=*(p_saved_amplitudes);
  }
}

Hadron_Decay_Handler* Hadron_Decays::ChooseDecayHandler(Particle* part)
{
#ifdef DEBUG__Hadrons
  if(part->Flav().IsStable()) {
    msg_Error()<<METHOD<<" for stable particle? part="<<part->Flav()<<endl;
  }
#endif
  for (HDHandlersIter hd=p_dechandlers->begin();hd!=p_dechandlers->end();hd++) {
    if (hd->second->CanDealWith(part->Flav().Kfcode())) {
      return hd->second;
    }
  }
  msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
    <<"   Unstable particle found ("<<part->Flav()<<") in "<<endl;
  if( part->ProductionBlob()->Type()==btp::Fragmentation ) {
    msg_Error()<<"   coming out of fragmentation."<<endl;
  }
  else {
    msg_Error()<<*(part->ProductionBlob())<<endl;
  }
  msg_Error()<<"   but no handler found to deal with it."<<std::endl
    <<"   Will continue and hope for the best."<<std::endl;
  return NULL;
}

void Hadron_Decays::CleanUp() {}

void Hadron_Decays::Finish(const std::string &) {
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
  map<kf::code,TH1D*>::iterator it;
  for(it=mass_hists.begin();it!=mass_hists.end();it++) {
    if(it->second->GetEntries()>0) {
      it->second->Write();
    }
  }
  p_file->Close();
#endif
#endif
}
