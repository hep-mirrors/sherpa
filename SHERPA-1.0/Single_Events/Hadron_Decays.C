#include "Hadron_Decays.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Mass_Handler.H"
#include "Momenta_Stretcher.H"

#include <utility>
#include <algorithm>

#ifdef PROFILE__Hadron_Decays
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

#ifdef DEBUG
#ifdef USING__Hadrons
#include "Hadron_Decay_Channel.H"
#endif
#include "TH1D.h"
#include "TFile.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadron_Decays::Hadron_Decays(HDHandlersMap * _dechandlers) :
  p_dechandlers(_dechandlers)
{
#ifdef DEBUG
  Fl_Iter fli;
  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.IsOn() && flav.IsHadron()) {
      double min = flav.PSMass()-3.0*flav.Width();
      double max = flav.PSMass()+3.0*flav.Width();
      p_file = new TFile("masses.root","RECREATE");
      TH1D myhist = TH1D(flav.ShellName().c_str(),flav.IDName().c_str(),100,min,max);
      myhist.SetDirectory(p_file);
      mass_hists.insert(make_pair(flav.Kfcode(),myhist));
    }
  }
#endif
  m_name      = std::string("Hadron_Decays");
  m_type      = eph::Hadronization;
}

Hadron_Decays::~Hadron_Decays()
{
}

bool SortByWidth(Particle* p1, Particle* p2) {
  return p1->Flav().Width() < p2->Flav().Width();
}

Return_Value::code Hadron_Decays::Treat(ATOOLS::Blob_List * bloblist, double & weight) 
{
  msg.Tracking()<<"--------- "<<METHOD<<" - START --------------"<<endl;
  PROFILE_HERE;
  if(p_dechandlers->empty()) return Return_Value::Nothing;

  if (bloblist->empty()) {
    msg.Error()<<"Potential error in "<<METHOD<<endl
      <<"   Incoming blob list contains "<<bloblist->size()<<" entries."<<endl
      <<"   Continue and hope for the best."<<endl;
    return Return_Value::Nothing;
  }

  bool found(true), didit(false);
  while (found) {
    found = false;
    Blob_List::iterator blit;
    for (blit=bloblist->begin();blit!=bloblist->end();++blit) {
      if ((*blit)->Has(blob_status::needs_hadrondecays)) {
        Particle_Vector daughters = (*blit)->GetOutParticles();
        // Sort daughters by width, necessary for mass treatment
        sort(daughters.begin(), daughters.end(), SortByWidth);

        // Preparing for mass smearing
        vector<double> widths(daughters.size());
        vector<double> masses(daughters.size());
        vector<Vec4D> saved_momenta; // needed for stretching/boosting back later
        vector<double> max_masses(daughters.size());
        for(size_t i=0;i<daughters.size();i++) {
          widths[i]=daughters[i]->Flav().Width();
          masses[i]=daughters[i]->Flav().PSMass(); // temporary storage of peak mass
          saved_momenta.push_back(daughters[i]->Momentum());
        }
        Vec4D total(0.0,0.0,0.0,0.0);
        for(int i=0;i<(*blit)->NInP();i++) total += (*blit)->InParticle(i)->Momentum();
        double max_mass = total.Mass();

        // Decaying all daughters one by one (in CMS*)
        bool retry_all = false;
        int all_again_trials = 0;
        do { // retry smearing of all daughters in case a single-mass retry didn't succeed
          retry_all = false;
          for (size_t i=0;i<daughters.size();i++) {
            max_masses[i] = i>0 ? max_masses[i-1]-masses[i-1] : max_mass;
            Particle * part = daughters[i];
            Mass_Handler masshandler(part->Flav());
            if( part->Status()==part_status::active && !part->Flav().IsStable() ) {
              Hadron_Decay_Handler * hdhandler = NULL;
              for (HDHandlersIter hd=p_dechandlers->begin();hd!=p_dechandlers->end();hd++) {
                if (hd->second->CanDealWith(part->Flav().Kfcode())) {
                  hdhandler = hd->second;
                  break;
                }
              }
              if (hdhandler) {
                Return_Value::code ret = Return_Value::Success;
                int trials=0;
                do { // retry decay with new m* as long as Retry_Method
                  trials++;
                  masses[i] = masshandler.GetMass(0,max_masses[i]);
                  part->SetFinalMass(masses[i]);
                  part->SetMomentum(Vec4D(masses[i],0.0,0.0,0.0));
                  Blob* blob;
                  // check if particle has a decayblob already, where its
                  // decay channel was stored in a previous try (HADRONS)
                  if(part->DecayBlob()) blob = part->DecayBlob();
                  else {
                    blob = new Blob();
                    blob->SetId();
                    blob->SetType(btp::Hadron_Decay);
                    blob->SetStatus(blob_status::needs_hadrondecays);
                    blob->AddToInParticles(part);
                    bloblist->push_back(blob);
                  }
                  ret = hdhandler->FillHadronDecayBlob(blob);
                  if( ret == Return_Value::Success ) found = true;
                  else if (ret == Return_Value::Retry_Method) {
                    if(trials>100) {
#ifdef DEBUG
                      msg.Error()<<endl<<"Warning in "<<METHOD<<" event "
                        <<rpa.gen.NumberOfDicedEvents()<<endl
                        <<"   Hadron_Decay_Handler "<<hdhandler->Name()<<" rejected mass of "<<endl
                        <<(*part)<<endl
                        <<"   retried method too often. Will retry all masses."<<endl;
                      Blob_Data_Base* data = (*blob)["hdc"];
                      if(data) {
                        HADRONS::Hadron_Decay_Channel *hdc;
                        hdc = (*blob)["hdc"]->Get<HADRONS::Hadron_Decay_Channel*>();
                        msg.Error()<<"   Chosen decay channel was "<<hdc->Name()<<endl;
                      }
                      msg.Error()<<"   Maximum masses so far:"<<endl;
                      for(size_t j=0;j<=i;j++) {
                        msg.Error()<<"   max_masses["<<j<<"]="<<max_masses[j]<<" for "
                          <<daughters[j]->Flav()<<" -> diced "<<masses[j]<<endl;
                      }
                      msg.Error()<<(**blit)<<endl; // show preliminary blob
#endif
                      for(size_t i=0;i<daughters.size();i++) {
                        daughters[i]->SetStatus(part_status::active);
                        bloblist->Delete(daughters[i]->DecayBlob());
                      }
                      retry_all = true;
                      all_again_trials++;
                      break;
                    }
                  }
                  else if (ret == Return_Value::Nothing) {
                    bloblist->Delete(blob);
                    part->SetStatus(part_status::active);
                  }
                  else {
                    msg.Error()<<"Error in "<<METHOD<<":"<<endl
                      <<"   Hadron_Decay_Handler "<<hdhandler->Name()<<" failed to decay "<<endl
                      <<"   "<<(*part)<<","<<endl
                      <<"   it returned "<<ret<<". Will retry event."<<endl;
                    return Return_Value::Retry_Event;
                  }
                } while( ret == Return_Value::Retry_Method );
              }
              else { // if no decay handler found (-> particle quasi stable)
                msg.Error()<<"Warning in "<<METHOD<<":"<<std::endl
                  <<"   Unstable particle found ("<<part->Flav()<<"), "<<endl;
                if( (*blit)->Type()==btp::Fragmentation ) {
                  msg.Error()<<"   coming out of fragmentation."<<endl;
                }
                else {
                  msg.Error()<<*(part->ProductionBlob())<<endl;
                }
                msg.Error()<<"   but no handler found to deal with it."<<std::endl
                  <<"   Will continue and hope for the best."<<std::endl;
                masses[i] = masshandler.GetMass(0,max_masses[i]);
                part->SetFinalMass(masses[i]);
                part->SetMomentum(Vec4D(masses[i],0.0,0.0,0.0)); // actually not needed
              }
            }
            else { // if particle stable
              masses[i] = masshandler.GetMass(0,max_masses[i]);
              part->SetFinalMass(masses[i]);
              part->SetMomentum(Vec4D(masses[i],0.0,0.0,0.0)); // actually not needed
            }
            if( retry_all ) break;
          }
          if( all_again_trials > 5 ) {
            msg.Error()<<"Warning in "<<METHOD<<endl
              <<"   retried mass dicing too often and didn't succeed. "
              <<"Will retry event."<<endl;
            return Return_Value::Retry_Event;
          }
        } while( retry_all == true ) ;
        // by here all daughters should be on their new mass shells and have either
        // a decayblob or momentum, in CMS.
        
        // now we stretch the initial momenta to the accepted m*
        Momenta_Stretcher stretch;
        stretch.StretchMomenta( daughters, saved_momenta);

        // now we boost all decayblobs back to these stretched target momenta
        for (size_t i=0;i<daughters.size();i++) {
          if(daughters[i]->DecayBlob()) {
            Poincare boost(saved_momenta[i]);
            boost.Invert();
            daughters[i]->DecayBlob()->Boost(boost);
            daughters[i]->SetStatus(part_status::decayed);
#ifdef DEBUG // fill the mass histogram of this flavour
            mass_hists[daughters[i]->Flav().Kfcode()].Fill(daughters[i]->FinalMass());
#endif
          }
          else {
            daughters[i]->SetMomentum(saved_momenta[i]);
#ifdef DEBUG // fill the mass histogram of this flavour
            mass_hists[daughters[i]->Flav().Kfcode()].Fill(daughters[i]->FinalMass());
#endif
          }
        }

        // finally set the position and lifetime in case it's a hadron decay blob
        if((*blit)->Type() == btp::Hadron_Decay) {
          double time        = (*blit)->InParticle(0)->LifeTime();
          Vec3D      spatial = (*blit)->InParticle(0)->Distance( time );
          Vec4D     position = Vec4D( time*rpa.c(), spatial );
          (*blit)->SetPosition( (*blit)->InParticle(0)->XProd() + position );
        }
        (*blit)->UnsetStatus(blob_status::needs_hadrondecays);

        didit = true;
      }
    }
  }
  msg.Tracking()<<"--------- Hadron_Decays::Treat - FINISH -------------"<<endl;
  return (didit ? Return_Value::Success : Return_Value::Nothing);
}

void Hadron_Decays::CleanUp() {}

void Hadron_Decays::Finish(const std::string &) {
#ifdef DEBUG
//   TFile myfile("masses.root","RECREATE");
  map<kf::code,TH1D>::iterator it;
  for(it=mass_hists.begin();it!=mass_hists.end();it++) {
    if(it->second.GetEntries()>0) {
      it->second.Write();
    }
  }
//   myfile.Close();
#endif
}
