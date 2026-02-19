#include "RECONNECTIONS/Main/Reconnection_Handler.H"
#include "RECONNECTIONS/Main/Reconnect_By_Singlet.H"
#include "RECONNECTIONS/Main/Reconnect_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Handler::Reconnection_Handler(const bool & on) :
  m_on(on),
  p_reconnector(new Reconnect_Statistical()),
  m_nfails(0)
{}

Reconnection_Handler::~Reconnection_Handler() {
  if (m_on) {
    msg_Info()<<METHOD<<": reconnection handler winds down with "
	      <<m_nfails<<" errors overall.\n";
  }
  delete p_reconnector;
}

void Reconnection_Handler::Initialize() {
  if (m_on) p_reconnector->Initialize();
}

void Reconnection_Handler::Reset() {
  if (m_on) p_reconnector->Reset();
}

Return_Value::code Reconnection_Handler::operator()(Blob_List *const blobs,
						    Particle_List *const parts) {
  if (!m_on) return Return_Value::Nothing;
  switch ((*p_reconnector)(blobs)) {
  case -1:
    // things went wrong, try new event and hope it works better
    if (m_nfails<5)
      msg_Error()<<"Error in "<<METHOD<<": reconnections didn't work out.\n"
		 <<"   Ask for new event and hope for the best.\n";
    p_reconnector->Reset();
    m_nfails++;
    return Return_Value::New_Event;
  case 1:
    // added colour reconnections, but produce a reconnection blob.
    AddReconnectionBlob(blobs);
  case 0:
    // didn't find any blob that needed reconnections
    break;
  }
  p_reconnector->Reset();

  auto variation_weights = p_reconnector->GetVariationWeights();
  
  if (auto* statistical_reconnector = dynamic_cast<Reconnect_Statistical*>(p_reconnector)) {
    const double weight_cutoff = statistical_reconnector->GetWeightCutoff();
    if (weight_cutoff > 0.) {
      for (size_t i = 0; i < variation_weights.size(); ++i) {
        if (variation_weights[i] > weight_cutoff) {
          variation_weights[i] = weight_cutoff;
        }
      }
    }
  }
  
  Blob *blob(blobs->FindFirst(btp::Signal_Process));
  if (blob == NULL)
    blob = blobs->FindFirst(btp::Hard_Collision);
  auto &wgt_map = (*blob)["WeightsMap"]->Get<Weights_Map>();

  const size_t n_rec_variations = variation_weights.size();

  auto it = wgt_map.find("MPI+RECONNECTIONS");
  const bool has_mpi = (it != wgt_map.end());
  const size_t n_mpi_variations = has_mpi ? (it->second.Size() - 1) : 0;

  const size_t n_max = std::max(n_mpi_variations, n_rec_variations);
  auto& wgt_category = wgt_map["MPI+RECONNECTIONS"];
  
  for (size_t i = 0; i < n_max; ++i) {
    const std::string name = "v" + std::to_string(i);
    const double mpi_wgt = (i < n_mpi_variations) ? wgt_category[name] : 1.;
    const double rec_wgt = (i < n_rec_variations) ? variation_weights[i] : 1.;
    wgt_category[name] = mpi_wgt * rec_wgt;
  }

  p_reconnector->ResetVariationWeights(n_rec_variations);

  return Return_Value::Success; 
}

void Reconnection_Handler::AddReconnectionBlob(Blob_List *const blobs) {
  Blob * blob = new Blob();
  blob->AddStatus(blob_status::needs_hadronization);
  blob->SetType(btp::Fragmentation);
  blob->SetId();
  Part_List * particles = p_reconnector->GetParticles();
  while (!particles->empty()) {
    Particle * part = particles->front();
    part->DecayBlob()->AddToOutParticles(part);
    part->SetDecayBlob(NULL);
    blob->AddToInParticles(part);
    particles->pop_front();
  }
  blobs->push_back(blob);
}

