#include "RECONNECTIONS/Main/Reconnection_Reweighting.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <algorithm>
#include <cmath>
#include "ATOOLS/Org/MyStrStream.H" // OUTPUT

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Reweighting::Reconnection_Reweighting() :
  m_n_variations(1),
  m_max_reweight_factor(-1.),

  // OUTPUT
  m_reweighting_output(false),
  p_blobs(NULL),
  m_total_events(0),
  m_n_reconnection_count(0),
  m_event_n_reconnections(0)
{}

Reconnection_Reweighting::~Reconnection_Reweighting() {
  // OUTPUT
  PrintVariationStatistics();
  if (m_cr_weight_file.is_open()) {
    m_cr_weight_file.close();
  }
  if (m_total_weight_file.is_open()) {
    m_total_weight_file.close();
  }
}

void Reconnection_Reweighting::Initialize() {
  auto s = Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
  m_etaQ2     = s["ETA_Q"].SetDefault({0.63}).GetVector<double>();
  m_reshuffle = s["RESHUFFLE"].SetDefault({1./9.}).GetVector<double>();

  for (size_t i{0}; i<m_etaQ2.size(); ++i)
    m_etaQ2[i] = sqr(m_etaQ2[i]);

  m_n_variations = std::max({m_etaQ2.size(), m_reshuffle.size(), size_t(1)});
  m_etaQ2.resize(m_n_variations, m_etaQ2[0]);
  m_reshuffle.resize(m_n_variations, m_reshuffle[0]);

  m_max_reweight_factor = s["MAX_REWEIGHT_FACTOR"].SetDefault(-1.).Get<double>();
  ResetEvent();

  // OUTPUT
  ResetStats();
  m_reweighting_output = s["REWEIGHTING_OUTPUT"].SetDefault(0).Get<int>() > 0;
  if (m_n_variations > 1) {
    m_cutoff_count.resize(m_n_variations, 0);
    m_sum_weights.resize(m_n_variations, 0.0);
    m_sum_weights_squared.resize(m_n_variations, 0.0);
  }
  m_total_events = 0;
  if (m_reweighting_output) {
    std::string filename = "cr_weights.dat";
    m_cr_weight_file.open(filename);
    std::string total_filename = "total_weights.dat";
    m_total_weight_file.open(total_filename);
  }
  if (m_cr_weight_file.is_open()) {
    m_cr_weight_file << "# n_cr";
    for (size_t ivar=1; ivar<m_n_variations; ivar++) {
      m_cr_weight_file << " w_v" << ivar;
    }
    m_cr_weight_file << "\n";
    m_cr_weight_file << std::scientific << std::setprecision(10);
  }
  if (m_total_weight_file.is_open()) {
      m_total_weight_file << "# w_total \n";
    m_total_weight_file << std::scientific << std::setprecision(10);
  }
}

void Reconnection_Reweighting::ResetEvent() {
  m_variation_weights.resize(m_n_variations);
  std::fill(m_variation_weights.begin(), m_variation_weights.end(), 1.);
}

void Reconnection_Reweighting::AcceptRejectReweighting(bool accepted,
                                const std::vector<double>& probs) {
  if (accepted) ++m_n_reconnection_count; // OUTPUT
  if (m_n_variations <= 1) return;
  if (probs[0] <= 1e-8) return;
  if (accepted) {
    // Accepted: weight *= p_var / p_nom
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      m_variation_weights[ivar] *= probs[ivar] / probs[0];
    }
  } else {
    // Rejected: weight *= (1 - p_var) / (1 - p_nom)
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      m_variation_weights[ivar] *= (1. - probs[ivar]) / (1. - probs[0]);
    }
  }
}

void Reconnection_Reweighting::ApplyVariationWeights(Blob_List *const blobs) {
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    double w_total = m_variation_weights[ivar];
    if (m_max_reweight_factor > 0. && w_total > m_max_reweight_factor) {
      w_total = m_max_reweight_factor;
      m_cutoff_count[ivar]++; // OUTPUT
    }
    m_variation_weights[ivar] = w_total;
  }
  Blob *blob(blobs->FindFirst(btp::Signal_Process));
  if (blob == NULL) blob = blobs->FindFirst(btp::Hard_Collision);
  if (blob != NULL) {
    auto wgtmap = (*blob)["WeightsMap"]->Get<Weights_Map>();
    CombineSoftPhysicsVariations(wgtmap, m_variation_weights);
    blob->AddData("WeightsMap", new Blob_Data<Weights_Map>(wgtmap));
    AccumulateEventStatistics(); // OUTPUT
  }
  ResetEvent();
}

///////////////////////////// OUTPUT AND STATISTICS /////////////////////////////

void Reconnection_Reweighting::CacheBlobs(Blob_List * blobs) {
  if (p_blobs == NULL) p_blobs = blobs;
}

void Reconnection_Reweighting::ResetStats() {
  p_blobs = NULL;
  m_n_reconnection_count = 0;
  m_event_variation_weights.resize(m_n_variations);
  std::fill(m_event_variation_weights.begin(),
            m_event_variation_weights.end(), 1.);
  m_event_n_reconnections = 0;
}

void Reconnection_Reweighting::AccumulateEventStatistics() {
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    m_event_variation_weights[ivar] *= m_variation_weights[ivar];
  }
  m_event_n_reconnections += m_n_reconnection_count;
  m_n_reconnection_count = 0;
}

void Reconnection_Reweighting::WriteEventStatistics() {
  if (m_cr_weight_file.is_open()) {
    m_cr_weight_file << m_event_n_reconnections;
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      m_cr_weight_file << " " << m_event_variation_weights[ivar];
    }
    m_cr_weight_file << "\n";
  }
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    const double w = m_event_variation_weights[ivar];
    m_sum_weights[ivar] += w;
    m_sum_weights_squared[ivar] += w * w;
  }
  Blob *blob(p_blobs->FindFirst(btp::Signal_Process));
  if (blob == NULL) blob = p_blobs->FindFirst(btp::Hard_Collision);
  if (blob != NULL) {
    auto wgtmap = (*blob)["WeightsMap"]->Get<Weights_Map>();
    if (m_total_weight_file.is_open()) {
      auto sw_it = wgtmap.find(SoftPhysicsKey);
      if (sw_it != wgtmap.end()) {
        const auto& sw = sw_it->second;
        for (size_t i=1; i<sw.Size(); ++i) {
          m_total_weight_file << " " << sw[i];
        }
      }
      m_total_weight_file << "\n";
    }
  }
  m_total_events++;
  ResetStats();
}

void Reconnection_Reweighting::PrintVariationStatistics() {
  if (m_n_variations <= 1 || m_total_events == 0) return;

  const std::string title = "CR Reweighting Statistics (events: " +
                            ToString<size_t>(m_total_events) + ") ";

  const int variation_col_size = std::max<int>(
      static_cast<int>(std::string("Variation").size()),
      static_cast<int>(("v" + ToString<size_t>(m_n_variations - 1)).size()));
  const int col_size = 15;

  int table_size = variation_col_size + col_size + col_size + col_size + 4;
  if (m_max_reweight_factor > 0.0) {
    table_size += col_size + col_size;
  }
  table_size = std::max(table_size, static_cast<int>(title.size()) + 4);

  msg_Out() << Frame_Header{table_size};
  MyStrStream line;
  line << om::bold << std::left << title << om::reset;
  msg_Out() << Frame_Line{line.str(), table_size};

  msg_Out() << Frame_Separator{table_size};
  line.str("");
  line << std::left << std::setw(variation_col_size) << "Variation"
       << std::right << std::setw(col_size) << "Avg. weight"
       << std::right << std::setw(col_size) << "ESS"
       << std::right << std::setw(col_size) << "ESS ratio";
  if (m_max_reweight_factor > 0.0) {
    line << std::right << std::setw(col_size) << "Cutoffs"
         << std::right << std::setw(col_size) << "Cutoff ratio";
  }
  msg_Out() << Frame_Line{line.str(), table_size};
  msg_Out() << Frame_Separator{table_size};

  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    const double sum_w = m_sum_weights[ivar];
    const double sum_w2 = m_sum_weights_squared[ivar];
    const double avg_weight = (m_total_events > 0) ? sum_w / m_total_events : 0.;
    const double ess = (sum_w2 > 0.) ? (sum_w * sum_w) / sum_w2 : 0.;
    const double ess_ratio = (m_total_events > 0) ? ess / m_total_events : 0.;
    const double cutoff_ratio =
        (m_total_events > 0) ? static_cast<double>(m_cutoff_count[ivar]) / m_total_events : 0.;

    line.str("");
    line << om::bold << std::left << std::setw(variation_col_size)
         << ("v" + ToString<size_t>(ivar)) << om::reset
         << std::right << om::brown << std::setw(col_size)
         << std::fixed << std::setprecision(6) << avg_weight
         << std::setw(col_size)
         << std::fixed << std::setprecision(1) << ess
         << std::setw(col_size)
         << std::fixed << std::setprecision(6) << ess_ratio << om::reset;
    if (m_max_reweight_factor > 0.0) {
      line << om::red << std::setw(col_size)
           << std::fixed << std::setprecision(1) << m_cutoff_count[ivar]
           << std::setw(col_size)
           << std::fixed << std::setprecision(6) << cutoff_ratio << om::reset;
    }
    msg_Out() << Frame_Line{line.str(), table_size};
  }
  msg_Out() << Frame_Footer{table_size};
}
