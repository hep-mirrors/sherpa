#include "ALPACA/Main/Alpaca.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


Alpaca::Alpaca(){   
    HIPars.Init();
    p_generator = std::make_shared<Event_Generator>();
    p_partons = std::make_shared<std::list<std::shared_ptr<Parton>>>();
    //m_translator(P2P_Translator());
}

bool Alpaca::operator()(ATOOLS::Blob_List * blobs) {
  // check if list isn't empty - we may have to delete them first.
  p_partons->clear();
  m_inparticles.clear();
  msg_Out() << "ALPACA called with blob" << endl;
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
          //msg_Out() << "Adding partons" << endl;
        }
      }
    } else{
      //msg_Out() << "Blob does NOT have blob_status::needs_rescattering" << endl;
    }
  } 
  if (p_partons->empty()){
    msg_Out() << "No partons added from blobs" << endl;
    return false;
  }
  for (std::list<std::shared_ptr<Parton>>::iterator piter = p_partons->begin();
    piter != p_partons->end(); ++piter) {
    //msg_Out() << **piter << std::endl;
}
  return true;
}

bool Alpaca::AddBlob(ATOOLS::Blob_List * blobs) {
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
    msg_Out() << "Adding parton " << *(p_partons->front()) << " to blob" << endl;
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
  msg_Out() << METHOD << endl;
    if(p_generator){
      msg_Out() << METHOD << " Evolve" << endl;
      p_generator->GenerateEvent(p_partons);
    } else{
      msg_Out() << METHOD << ": ERROR: parton list not initialized, will exit" << endl;
      exit(1.);
    }
    return true;
}
