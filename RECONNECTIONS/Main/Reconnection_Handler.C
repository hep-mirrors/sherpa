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
  m_nfails(0),
  m_n_variations(0),
  m_n_mpi_variations(SIZE_MAX),
  m_total_events(0)
{}

Reconnection_Handler::~Reconnection_Handler() {
  if (m_on) {
    msg_Info()<<METHOD<<": reconnection handler winds down with "
	      <<m_nfails<<" errors overall.\n";
    PrintVariationStatistics();
  }
  // output
  if (m_cr_weight_file.is_open()) {
    m_cr_weight_file.close();
  }
  if (m_total_weight_file.is_open()) {
    m_total_weight_file.close();
  }
  // output
  delete p_reconnector;
}

void Reconnection_Handler::Initialize() {
  if (m_on) p_reconnector->Initialize();

  m_n_variations = p_reconnector->GetVariationSize();
  m_weight_cutoff = p_reconnector->GetWeightCutoff();
  m_cutoff_count.resize(m_n_variations, 0);
  m_sum_weights.resize(m_n_variations, 0.0);
  m_sum_weights_squared.resize(m_n_variations, 0.0);

  p_reconnector->ResetVariationWeights(m_n_variations);

  // output
  const bool print_files = p_reconnector->GetWeightOutput();
  if (print_files && m_on) {
    std::string filename = "cr_weights.dat";
    m_cr_weight_file.open(filename);
    std::string total_filename = "total_weights.dat";
    m_total_weight_file.open(total_filename);
  }
  
  if (m_cr_weight_file.is_open()) {
    m_cr_weight_file << "# n_cr";
    for (size_t i=0; i<m_n_variations; i++) {
      m_cr_weight_file << " w_cr_v" << i;
    }
    m_cr_weight_file << "\n";
    m_cr_weight_file << std::scientific << std::setprecision(10);
  }
  if (m_total_weight_file.is_open()) {
      m_total_weight_file << "#";
      m_total_weight_file << " w_total";
      m_total_weight_file << "\n";
    m_total_weight_file << std::scientific << std::setprecision(10);
  }
  // output
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
    p_reconnector->ResetVariationWeights(m_n_variations);
    m_nfails++;
    return Return_Value::New_Event;
  case 1:
    // added colour reconnections, but produce a reconnection blob.
    AddReconnectionBlob(blobs);
  case 0:
    // didn't find any blob that needed reconnections
    break;
  }

  auto variation_weights = p_reconnector->GetVariationWeights();
  size_t n_reconnections = p_reconnector->GetReconnectionCount();

  msg_Out() << METHOD << ": " << n_reconnections << " WEIHGTS: ";
  for (size_t i = 0; i < m_n_variations; ++i) {
    msg_Out() << variation_weights[i] << " ";
  }
  msg_Out() << "\n";
  
  p_reconnector->Reset();
  
  if (m_weight_cutoff > 0.) {
    for (size_t i = 0; i < m_n_variations; ++i) {
      if (variation_weights[i] > m_weight_cutoff) {
        variation_weights[i] = m_weight_cutoff;
        m_cutoff_count[i]++;
      }
    }
  }
  
  Blob *blob(blobs->FindFirst(btp::Signal_Process));
  if (blob == NULL)
    blob = blobs->FindFirst(btp::Hard_Collision);
  auto &wgt_map = (*blob)["WeightsMap"]->Get<Weights_Map>();

  auto it = wgt_map.find("MPI+RECONNECTIONS");
  
  size_t n_mpi_variations;
  if (m_n_mpi_variations == SIZE_MAX) {
    const bool has_mpi = (it != wgt_map.end());
    n_mpi_variations = has_mpi ? (it->second.Size() - 1) : 0;
    m_n_mpi_variations = n_mpi_variations;
  } else {
    n_mpi_variations = m_n_mpi_variations;
  }

  const size_t n_max = std::max(n_mpi_variations, m_n_variations);
  auto& wgt_category = wgt_map["MPI+RECONNECTIONS"];
  
  std::vector<double> combined_weights(n_max);
  for (size_t i = 0; i < n_max; ++i) {
    const std::string name = "v" + std::to_string(i);
    const double w_mpi = (i < n_mpi_variations) ? wgt_category[name] : 1.;
    const double w_rec = (i < m_n_variations) ? variation_weights[i] : 1.;
    combined_weights[i] = w_mpi * w_rec;
    wgt_category[name] = combined_weights[i];
  }
  
  // output
  if (m_cr_weight_file.is_open()) {
    m_cr_weight_file << n_reconnections;
    for (size_t i = 0; i < m_n_variations; ++i) {
      m_cr_weight_file << " " << variation_weights[i];
    }
    m_cr_weight_file << "\n";
  // output
  }
  if (m_total_weight_file.is_open()) {
    for (size_t i = 0; i < n_max; ++i) {
      if (i > 0) m_total_weight_file << " ";
      m_total_weight_file << combined_weights[i];
    }
    m_total_weight_file << "\n";
  }
  
  m_total_events++;
  for (size_t i = 0; i < m_n_variations; ++i) {
    double w = variation_weights[i];
    m_sum_weights[i] += w;
    m_sum_weights_squared[i] += w * w;
  }

  p_reconnector->ResetVariationWeights(m_n_variations);

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

void Reconnection_Handler::PrintVariationStatistics() {
  if (m_n_variations > 1 && m_total_events > 0) {
    msg_Info() << "\n" 
               << "   " << std::string(77, '-') << "\n"
               << "   | CR Reweighting Statistics: Total events: " 
               << std::setw(15) << m_total_events 
               << std::string(18,' ') << "|\n";
    
    msg_Info() << "   " << std::string(77, '-') << "\n"
               << "   | " 
               << std::setw(8) << "v#" 
               << " | " << std::setw(15) << "ESS"
               << " | " << std::setw(10) << "ESS ratio"
               << " | " << std::setw(10) << "Cutoffs"
               << " | " << std::setw(18) << "Cutoffs ratio"
               << std::string(1,' ') << "|\n";

    msg_Info() << "   " << std::string(77, '-') << "\n";
    for (size_t ivar = 0; ivar < m_n_variations; ++ivar) {
      double sum_w = m_sum_weights[ivar];
      double sum_w2 = m_sum_weights_squared[ivar];
      double ess = (sum_w2 > 0) ? (sum_w * sum_w) / sum_w2 : 0.0;
      double ess_ratio = (m_total_events > 0) ? ess / m_total_events : 0.0;
      double cutoff_ratio = (m_total_events > 0) ? static_cast<double>(m_cutoff_count[ivar]) / m_total_events : 0.0;
      
      msg_Info() << "   | " << std::setw(8) << ivar
                 << " | " << std::setw(15) << std::fixed << std::setprecision(1) << ess
                 << " | " << std::setw(10) << std::fixed << std::setprecision(6) << ess_ratio
                 << " | " << std::setw(10) << m_cutoff_count[ivar]
                 << " | " << std::setw(10) << std::fixed << std::setprecision(6) << cutoff_ratio
                 << std::string(9,' ') << "|\n";
    }
    msg_Info() << "   " << std::string(77, '-') << "\n";
  }
}
