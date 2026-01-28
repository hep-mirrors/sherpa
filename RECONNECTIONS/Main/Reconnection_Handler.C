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

  const auto& variation_weights = p_reconnector->GetVariationWeights();
  Blob *blob(blobs->FindFirst(btp::Signal_Process));
  if (blob == NULL)
    blob = blobs->FindFirst(btp::Hard_Collision);
  auto &wgt_map = (*blob)["WeightsMap"]->Get<Weights_Map>();

  const size_t n_variations = variation_weights.size();
  const size_t n_mpi_variations = wgt_map["MPI"].Size() - 1;

  // static constexpr double max_weight = 2e2;
  // for (size_t i=0; i<n_variations; i++) {
  //     variation_weights[i] = std::min(variation_weights[i], max_weight);
  // }

  wgt_map["RECONNECTIONS"]["v0"] = variation_weights[0];
  if (n_mpi_variations>0) wgt_map["MPI+RECONNECTIONS"]["v0"] = wgt_map["MPI"]["v0"] * variation_weights[0];

  const size_t n_combined = std::min(n_mpi_variations, n_variations);
  for (size_t i=1; i<n_combined; ++i) {
    const std::string name = "v" + std::to_string(i);
    wgt_map["MPI+RECONNECTIONS"][name] = wgt_map["MPI"][name] * variation_weights[i];
  }
  for (size_t i=n_combined==0 ? 1 : n_combined; i<n_variations; ++i) {
    const std::string name = "v" + std::to_string(i);
    wgt_map["RECONNECTIONS"][name] = variation_weights[i];
  }

  if (n_mpi_variations>0) {
    const double tmp_mpi_nominal = wgt_map["MPI"]["v0"];
    std::vector<double> tmp_mpi_weights;
    if(n_mpi_variations > n_variations) {
      tmp_mpi_weights.resize(n_mpi_variations - n_variations, 1.0);
      for (size_t i=n_variations; i<n_mpi_variations; ++i) {
        const std::string name = "v" + std::to_string(i);
        tmp_mpi_weights[i - n_variations] = wgt_map["MPI"][name];
      }
    }
    wgt_map.erase("MPI");
    wgt_map["MPI"]["v0"] = tmp_mpi_nominal;
    if (n_mpi_variations > n_variations) {
      for (size_t i=n_variations; i<n_mpi_variations; ++i) {
        const std::string name = "v" + std::to_string(i);
        wgt_map["MPI"][name] = tmp_mpi_weights[i - n_variations];
      }
    }
  }

  p_reconnector->ResetVariationWeights(n_variations);

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

