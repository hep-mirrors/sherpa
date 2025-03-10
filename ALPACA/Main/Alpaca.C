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
}

bool Alpaca::Harvest(ATOOLS::Blob_List * blobs) {
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(blob_status::needs_rescattering)) {
      Particle_Vector partvec((*bit)->GetOutParticles());
      msg_Out()<<(**bit)<<"\n";
      for (vector<Particle *>::iterator piter=partvec.begin(); piter != partvec.end(); piter++) {
        m_partons.push_back(m_translator(*piter, (*bit)->Position()));
      }
    }
  }
  for (vector<Parton>::iterator piter=m_partons.begin(); piter != m_partons.end(); piter++) {
      msg_Out()<<*piter<<endl;
  }
  return false;
}
