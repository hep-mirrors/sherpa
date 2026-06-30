#include "RECONNECTIONS/Main/Reconnection_Reweighting.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <algorithm>
#include <cmath>

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Reweighting::Reconnection_Reweighting() :
  m_n_variations(1),
  m_max_reweight_factor(-1.)
{}

Reconnection_Reweighting::~Reconnection_Reweighting() {}

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
}

void Reconnection_Reweighting::ResetEvent() {
  m_variation_weights.resize(m_n_variations);
  std::fill(m_variation_weights.begin(), m_variation_weights.end(), 1.);
}

void Reconnection_Reweighting::AcceptRejectReweighting(bool accepted,
                                const std::vector<double>& probs) {
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
    }
    m_variation_weights[ivar] = w_total;
  }
  Blob *blob(blobs->FindFirst(btp::Signal_Process));
  if (blob == NULL) blob = blobs->FindFirst(btp::Hard_Collision);
  if (blob != NULL) {
    auto wgtmap = (*blob)["WeightsMap"]->Get<Weights_Map>();
    CombineSoftPhysicsVariations(wgtmap, m_variation_weights);
    blob->AddData("WeightsMap", new Blob_Data<Weights_Map>(wgtmap));
  }
  ResetEvent();
}
