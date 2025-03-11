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
      for (vector<Particle *>::iterator piter=outparts.begin(); piter != outparts.end(); piter++) {
        if((*piter)->Info() == 'F' ) m_partons.push_back(m_translator(*piter, (*bit)->Position()));
      }
    }
  }
  for (vector<Parton *>::iterator piter=m_partons.begin(); piter != m_partons.end(); piter++) {
      msg_Out()<<**piter<<endl;
  }
  return true;
}

bool Alpaca::AddBlob(ATOOLS::Blob_List * blobs) {
  Blob * blob = blobs->AddBlob(btp::Hard_Collision);
  blob->SetId();
  blob->SetTypeSpec("Rescattering");
  blob->SetStatus(blob_status::needs_reconnections);
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(blob_status::needs_rescattering)) {
      (*bit)->UnsetStatus(blob_status::needs_rescattering);
      vector<Particle *> outparts((*bit)->GetOutParticles());
      for (vector<Particle *>::iterator piter=outparts.begin(); piter != outparts.end(); piter++) {
        if((*piter)->Info() == 'F' ) blob->AddToInParticles(*piter);
      }
    }
  }
  for (vector<Parton *>::iterator piter = m_partons.begin(); piter != m_partons.end(); piter++) {
      blob->AddToOutParticles(m_translator((*piter)));
  }
  PRINT_VAR(*blobs);
  return true;
}

bool Alpaca::Evolve() {
    return true;
}
