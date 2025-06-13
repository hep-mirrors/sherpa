#include "ALPACA/Main/Alpaca.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "MODEL/Main/Model_Base.H"

#include <iostream>

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


Alpaca::Alpaca()
    //:m_mcekrt("runcards/mcekrt_001") //MC-EKRT, temporarily commented out
{   
    //msg_Out() << "FILENAME: runcards/mcekrt_" << to_string(ran->GetSeed()) << endl;
    HIPars.Init();
    p_partons = std::make_shared<std::list<std::shared_ptr<Parton>>>();
    p_flow_backtrack = std::make_shared<std::list<std::pair<int, int>>>();
    p_generator = std::make_shared<Event_Generator>(p_partons, p_flow_backtrack);
    //m_translator(P2P_Translator());
}

bool Alpaca::operator()(ATOOLS::Blob_List * blobs) {
  // check if list isn't empty - we may have to delete them first.
  p_partons->clear();
  p_flow_backtrack->clear();
  m_inparticles.clear();
  msg_Out() << "ALPACA called with blob" << endl;
  /* //MC-EKRT, temporarily commented out
  shared_ptr<list<shared_ptr<Parton>>> temp_list = m_mcekrt.GetPartons(1);
  msg_Out() << METHOD << ": Partons collected in ALPACA main from MC-EKRT" << endl;
  for (list<shared_ptr<Parton>>::iterator piter = temp_list->begin(); piter != temp_list->end(); ++piter) {
    msg_Out() << METHOD << **piter << std::endl;
  }
  */
  return (Harvest(blobs) && Evolve() && AddBlob(blobs));
}

bool Alpaca::Harvest(ATOOLS::Blob_List * blobs) {
  bool found_rescatter_blob = false;
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(blob_status::needs_rescattering)) {
      //msg_Out()<<(**bit)<<"\n";
      (*bit)->UnsetStatus(blob_status::needs_rescattering|
			  blob_status::needs_reconnections|
			  blob_status::needs_hadronization);
      vector<Particle *> outparts((*bit)->GetOutParticles());
      for (vector<Particle *>::iterator piter=outparts.begin(); piter != outparts.end(); piter++) {
        if ((*piter)->Status() == part_status::active &&
            (*piter)->Info() == 'F' &&
            (*piter)->Flav().Strong() &&
            !(*piter)->DecayBlob()) {
          m_inparticles.push_back((*piter));
          p_partons->push_back(m_translator(*piter));
        }
      }
    } else{
      //msg_Out() << "Blob does NOT have blob_status::needs_rescattering" << endl;
    }
  } 
  if (p_partons->empty()){
    msg_Out() << "WARNING: No partons added from blobs, no blob with blob_status::needs_rescattering in blob list" << endl;
    return false;
  } else{
    CheckIfSinglet();
    msg_Out() << METHOD << ": Partons have been added from blob, sending to ALPACA." << endl;
  }
  /*
  for (list<shared_ptr<Parton>>::iterator piter = p_partons->begin();
    piter != p_partons->end(); ++piter) {
    //msg_Out() << **piter << std::endl;
  }
  */
  return true;
}

bool Alpaca::AddBlob(ATOOLS::Blob_List * blobs) {
  //First backtrack colour from the evolution in ALPACA
  CheckFlowBacktrack(blobs);

  Blob * blob = blobs->AddBlob(btp::Hard_Collision);
  blob->SetId();
  blob->SetTypeSpec("ALPACA");
  blob->SetStatus(blob_status::needs_reconnections|
		  blob_status::needs_hadronization);
  for (vector<Particle *>::iterator piter = m_inparticles.begin();
       piter!=m_inparticles.end();piter++) {
    blob->AddToInParticles(*piter);
    (*piter)->SetStatus(part_status::decayed);
  }
  while (!p_partons->empty()) {
    //msg_Out() << "Adding parton " << *(p_partons->front()) << " to blob" << endl;
    blob->AddToOutParticles(m_translator(p_partons->front()));
    p_partons->pop_front(); 
  }
  //for (list<Parton *>::iterator piter = m_partons.begin();
  //   piter != m_partons.end(); piter++) {
  //  blob->AddToOutParticles(m_translator((*piter)));
  //}
  PRINT_VAR(*blobs);
  return true;
}

bool Alpaca::Evolve() {
    if(p_generator){
      msg_Out() << METHOD << " Evolve converted blob in ALPACA" << endl;
      p_generator->GenerateEvent();
      msg_Out() << METHOD << " Evolution in ALPACA complete, returning parton list" << endl;
    } else{
      msg_Out() << METHOD << ": ERROR: parton list not initialized, will exit" << endl;
      exit(1.);
    }
    return true;
}

void Alpaca::CheckIfSinglet(){
  int flow_id_1, flow_id_2;
  int found_id_1, found_id_2;
  int found_err_1, found_err_2;
  bool all_found = true;
  for (list<shared_ptr<Parton>>::iterator iter1=p_partons->begin(); iter1!=p_partons->end(); iter1++) {
    flow_id_1 = (*iter1)->GetFlow(1);
    flow_id_2 = (*iter1)->GetFlow(2);
    found_id_1 = 0;
    found_id_2 = 0;
    found_err_1 = 0;
    found_err_2 = 0;
    for (list<shared_ptr<Parton>>::iterator iter2=p_partons->begin(); iter2!=p_partons->end(); iter2++) {
      if((*iter1) == (*iter2)) continue;

      if(flow_id_1 == (*iter2)->GetFlow(2) && flow_id_1 != 0) found_id_1++;
      if(flow_id_2 == (*iter2)->GetFlow(1) && flow_id_2 != 0) found_id_2++;
      if(flow_id_1 == (*iter2)->GetFlow(1) && flow_id_1 != 0) found_err_1++;
      if(flow_id_2 == (*iter2)->GetFlow(2) && flow_id_2 != 0) found_err_2++;
    }

    if(found_id_1 != 1 && flow_id_1 != 0){
      msg_Out() << "WARNING: flow_id = " << flow_id_1 << " has " << found_id_1 << " anti-partners, adding to list" << endl;
      p_flow_backtrack->push_back(make_pair(flow_id_1, -1));
      all_found = false;
    }
    if(found_id_2 != 1 && flow_id_2 != 0){
      p_flow_backtrack->push_back(make_pair(flow_id_2, -1));
      msg_Out() << "WARNING: flow_id = " << flow_id_2 << " has " << found_id_2 << " anti-partners, adding to list" << endl;
      all_found = false;
    }
    if(found_err_1 != 0 && flow_id_1 != 0){
      msg_Out() << "WARNING: flow_id = " << flow_id_1 << " has " << found_err_1 << " copies" << endl;
      all_found = false;
    }
    if(found_err_2 != 0 && flow_id_2 != 0){
      msg_Out() << "WARNING: flow_id = " << flow_id_2 << " has " << found_err_2 << " copies" << endl;
      all_found = false;
    }
  }

  if(all_found){
    msg_Out() << "All colour indices match, total is singlet" << endl;
  } else{
    msg_Out() << "Error in colour indices, not a singlet in total" << endl;
    //exit(1.);
  }
}

void Alpaca::CheckFlowBacktrack(ATOOLS::Blob_List * blobs){
  msg_Out() << "Flow backtrack list" << endl;
  for (list<pair<int,int>>::iterator iter=p_flow_backtrack->begin(); iter!=p_flow_backtrack->end(); iter++) {
    msg_Out() << iter->first << " " << iter->second << endl;
    if(iter->second != -1){
      //Colour index have changed and needs to be backtracked in another blob
      bool found_index = false;
      unsigned int old_index = iter->first;
      unsigned int new_index = iter->second;
      for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
        vector<Particle *> outparts((*bit)->GetOutParticles());
        for (vector<Particle *>::iterator piter=outparts.begin(); piter != outparts.end(); piter++) {
          if ((*piter)->GetFlow(1) == old_index) {
            msg_Out() << "Backtracked external flow id = " << (*piter)->GetFlow(1) << " to " << new_index << endl;
            (*piter)->SetFlow(1, new_index);
            found_index = true;
          } else if ((*piter)->GetFlow(2) == old_index) {
            msg_Out() << "Backtracked external flow id = " << (*piter)->GetFlow(2) << " to " << new_index << endl;
            (*piter)->SetFlow(2, new_index);
            found_index = true;
          }

          if(found_index){
            break;
          }
        }

        if(found_index){
          break;
        }
      }

      if(!found_index){
        msg_Out() << METHOD << ": ERROR: trying to backtrack colour flow in blob external to ALPACA, cannot find index. Will exit." << endl;
        exit(1);
      }
    }
  }

  
}