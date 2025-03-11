#include "ALPACA/Main/Alpaca.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


Alpaca::Alpaca() {
    //m_translator(P2P_Translator());
}

bool Alpaca::operator()(ATOOLS::Blob_List * blobs) {
  // check if list isn't empty - we may have to delete them first.
  m_partons.clear();
  m_inparticles.clear();
  return (Harvest(blobs) && Evolve() && AddBlob(blobs));
}

bool Alpaca::Harvest(ATOOLS::Blob_List * blobs) {
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(blob_status::needs_rescattering)) {
      msg_Out()<<(**bit)<<"\n";
      (*bit)->UnsetStatus(blob_status::needs_rescattering|
			  blob_status::needs_reconnections|
			  blob_status::needs_hadronization);
      vector<Particle *> outparts((*bit)->GetOutParticles());
      for (vector<Particle *>::iterator piter=outparts.begin();
	   piter != outparts.end(); piter++) {
        if ((*piter)->Status() == part_status::active &&
	    (*piter)->Info() == 'F' &&
	    (*piter)->Flav().Strong() &&
	    !(*piter)->DecayBlob()) {
	  m_inparticles.push_back((*piter));
	  m_partons.push_back(m_translator(*piter));
	}
      }
    }
  }
  if (m_partons.empty()) return false;
  for (list<Parton *>::iterator piter=m_partons.begin();
       piter != m_partons.end(); piter++) {
      msg_Out()<<**piter<<endl;
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
  while (!m_partons.empty()) {
    blob->AddToOutParticles(m_translator(m_partons.front()));
    m_partons.pop_front();
  }
  //for (list<Parton *>::iterator piter = m_partons.begin();
  //   piter != m_partons.end(); piter++) {
  //  blob->AddToOutParticles(m_translator((*piter)));
  //}
  PRINT_VAR(*blobs);
  return true;
}

bool Alpaca::Evolve() {
    return true;
}
