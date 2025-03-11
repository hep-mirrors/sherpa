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
  if (!Harvest(blobs)) exit(1);
  if (!Evolve()) exit(1);
  if (!AddBlob(blobs)) exit(1);
  return true;
}

bool Alpaca::Harvest(ATOOLS::Blob_List * blobs) {
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(blob_status::needs_rescattering)) {
      msg_Out()<<(**bit)<<"\n";
      vector<Particle *> outparts((*bit)->GetOutParticles());
      for (vector<Particle *>::iterator piter=outparts.begin();
	   piter != outparts.end(); piter++) {
        if((*piter)->Status() == part_status::active &&
           (*piter)->Info() == 'F' &&
           (*piter)->Flav().Strong() &&
           !(*piter)->DecayBlob())
	      m_partons.push_back(m_translator(*piter));
      }
    }
  }
  for (vector<Parton *>::iterator piter=m_partons.begin();
       piter != m_partons.end(); piter++) {
      msg_Out()<<**piter<<endl;
  }
  return true;
}

bool Alpaca::AddBlob(ATOOLS::Blob_List * blobs) {
  Blob * blob = blobs->AddBlob(btp::Hard_Collision);
  blob->SetId();
  blob->SetTypeSpec("ALPACA");
  blob->SetStatus(blob_status::needs_reconnections|blob_status::needs_hadronization);
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(blob_status::needs_rescattering)) {
      (*bit)->UnsetStatus(blob_status::needs_rescattering|blob_status::needs_hadronization);
      vector<Particle *> outparts((*bit)->GetOutParticles());
      for (vector<Particle *>::iterator
	     piter=outparts.begin(); piter != outparts.end(); piter++) {
        if((*piter)->Info() == 'F' ) {
	  blob->AddToInParticles(*piter);
	  (*piter)->SetStatus(part_status::decayed);
	}
      }
    }
  }
  for (vector<Parton *>::iterator piter = m_partons.begin();
       piter != m_partons.end(); piter++) {
      blob->AddToOutParticles(m_translator((*piter)));
  }
  PRINT_VAR(*blobs);
  return true;
}

bool Alpaca::Evolve() {
    return true;
}
