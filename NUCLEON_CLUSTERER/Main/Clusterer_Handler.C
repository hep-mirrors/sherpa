#include "NUCLEON_CLUSTERER/Main/Clusterer_Handler.H"
#include "NUCLEON_CLUSTERER/Main/Cluster_By_Singlet.H"
#include "NUCLEON_CLUSTERER/Main/Cluster_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace NUCLEON_CLUSTERER;
using namespace ATOOLS;
using namespace std;

Clusterer_Handler::Clusterer_Handler(const bool & on) :
  m_on(on),
  p_clusterer(new Cluster_Statistical()),
  m_nfails(0)
{}

Clusterer_Handler::~Clusterer_Handler() {
  if (m_on) {
    msg_Info()<<METHOD<<": Clusterer handler winds down with "
	      <<m_nfails<<" errors overall.\n";
  }
  delete p_clusterer;
}

void Clusterer_Handler::Initialize() {
  if (m_on) p_clusterer->Initialize();
}

void Clusterer_Handler::Reset() {
  if (m_on) p_clusterer->Reset();
}

Return_Value::code Clusterer_Handler::operator()(Blob_List *const blobs,
						    Particle_List *const parts) {
  if (!m_on) return Return_Value::Nothing;
  switch ((*p_clusterer)(blobs)) {
  case -1:
    // things went wrong, try new event and hope it works better
    msg_Tracking()<<"Error in "<<METHOD<<": NUCLEON_CLUSTERER didn't work out.\n"
		  <<"   Ask for new event and hope for the best.\n";
    p_clusterer->Reset();
    m_nfails++;
    return Return_Value::New_Event;
  case 1:
    // added colour NUCLEON_CLUSTERER, but produce a Clusterer blob.
    AddClustererBlob(blobs);
  case 0:
    // didn't find any blob that needed NUCLEON_CLUSTERER
    break;
  }
  p_clusterer->Reset();
  return Return_Value::Success; 
}

void Clusterer_Handler::AddClustererBlob(Blob_List *const blobs) {
  Blob * blob = new Blob();
  blob->AddStatus(blob_status::needs_hadronization);
  blob->SetType(btp::Fragmentation);
  blob->SetId();
  Part_List * particles = p_clusterer->GetParticles();
  while (!particles->empty()) {
    Particle * part = particles->front();
    part->DecayBlob()->AddToOutParticles(part);
    part->SetDecayBlob(NULL);
    blob->AddToInParticles(part);
    particles->pop_front();
  }
  blobs->push_back(blob);
}

