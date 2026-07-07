#include "AHADIC++/Tools/Ahadic_Reweighting.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <algorithm>
#include <cmath>

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Ahadic_Reweighting::Ahadic_Reweighting() :
  m_n_variations(1),
  m_max_reweight_factor(-1.),
  m_reweight_max_nsplit(-1)
{}

Ahadic_Reweighting::~Ahadic_Reweighting() {}

void Ahadic_Reweighting::Initialize() {
  ///////////////////////////////////////////////////////////////////////////
  // Initialize the parameter variations for hadronization reweighting.
  ///////////////////////////////////////////////////////////////////////////
  static const std::vector<std::string> variation_keys = {
    "kT_0",
    "alphaG",
    "alphaL","betaL","gammaL",
    "alphaD","betaD","gammaD",
    "alphaB","betaB","gammaB",
    "alphaH","betaH","gammaH",
    "Strange_fraction","Baryon_fraction",
    "P_qs_by_P_qq","P_ss_by_P_qq","P_di_1_by_P_di_0"
  };
  m_n_variations = 1;
  for (const auto& key : variation_keys) {
    m_variation_vectors[key] = hadpars->GetVariationVector(key);
    m_n_variations = std::max(m_n_variations, m_variation_vectors[key].size());
  }
  for (const auto& key : variation_keys) {
    auto& vec = m_variation_vectors[key];
    vec.resize(m_n_variations, vec[0]);
  }
  for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
    const double strange = m_variation_vectors["Strange_fraction"][ivar];
    m_variation_vectors["P_qs_by_P_qq"][ivar] *= strange;
    m_variation_vectors["P_ss_by_P_qq"][ivar] *= sqr(strange);
  }

  m_max_reweight_factor = hadpars->Get("max_reweight_factor");
  m_reweight_max_nsplit = hadpars->Switch("reweight_max_nsplit");
  ResetEvent();
}

std::vector<double>
Ahadic_Reweighting::GetVariationVector(const std::string& keyword) const {
  ///////////////////////////////////////////////////////////////////////////
  // Access to the padded parameter-variation vectors for the tools that
  // compute the per-variation selection probabilities.
  ///////////////////////////////////////////////////////////////////////////
  auto piter = m_variation_vectors.find(keyword);
  if (piter != m_variation_vectors.end()) return piter->second;
  THROW(fatal_error, "Keyword " + keyword +
                     " not found in Ahadic_Reweighting vector map.");
}

void Ahadic_Reweighting::ResetEvent() {
  ///////////////////////////////////////////////////////////////////////////
  // Reset all variation weights to 1.0 for a new event.
  ///////////////////////////////////////////////////////////////////////////
  m_variation_weights.assign(m_n_variations, 1.);
  m_flavour_weights.assign(m_n_variations, 1.);
  m_gluon_weights.assign(m_n_variations, 1.);
  m_cluster_weights.assign(m_n_variations, 1.);
  m_soft_weights.assign(m_n_variations, 1.);
  m_kt_weights.assign(m_n_variations, 1.);
  m_tmp_flavour_weights.assign(m_n_variations, 1.);
  m_tmp_gluon_weights.assign(m_n_variations, 1.);
  m_tmp_cluster_weights.assign(m_n_variations, 1.);
  m_tmp_kt_weights.assign(m_n_variations, 1.);
}

void Ahadic_Reweighting::AcceptTmpWeights(std::vector<double>& weights,
                                       std::vector<double>& tmp_weights) {
  ///////////////////////////////////////////////////////////////////////////
  // Commit the weights of a successful splitting into the event weights.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    weights[ivar] *= tmp_weights[ivar];
  }
  std::fill(tmp_weights.begin(), tmp_weights.end(), 1.);
}

void Ahadic_Reweighting::ResetTmpWeights(std::vector<double>& tmp_weights) {
  ///////////////////////////////////////////////////////////////////////////
  // Discard the weights of a failed splitting attempt.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  std::fill(tmp_weights.begin(), tmp_weights.end(), 1.);
}

void Ahadic_Reweighting::FlavourSelectionReweighting(
    const std::vector<double>& popweights, const std::vector<double>& norms) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Flavour_Selector for each popped flavour.
  // The selection probability is p = popweight / norm, with the norm
  // depending on the event through the accessible flavours;
  // weight *= p_var / p_nom for the selected flavour.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  const double prob_nom = popweights[0] / norms[0];
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    m_tmp_flavour_weights[ivar] *= (popweights[ivar] / norms[ivar]) / prob_nom;
  }
}

void Ahadic_Reweighting::GluonSplittingReweighting(
    const std::vector<double>& probs) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Gluon_Splitter for each accepted z. The probs are the
  // normalised fragmentation-function densities f(z)/int f, so that the
  // accept/reject trials integrate out and rejected z play no role;
  // weight *= p_var / p_nom.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    m_tmp_gluon_weights[ivar] *= probs[ivar] / probs[0];
  }
}

void Ahadic_Reweighting::ClusterSplittingReweighting(
    const bool accepted, const std::vector<double>& probs) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Cluster_Splitter for each accept/reject trial of the z
  // selection, with each variation normalised to its own tight maximum;
  // probs[0] is the nominal acceptance probability.
  // For accepted z: weight *= p_var / p_nom
  // For rejected z: weight *= (1 - p_var) / (1 - p_nom)
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    if (accepted) {
      const double ratio = probs[ivar] / probs[0];
      if (!std::isnan(ratio)) {
        m_tmp_cluster_weights[ivar] *= ratio;
      }
    } else {
      m_tmp_cluster_weights[ivar] *= (1. - probs[ivar]) / (1. - probs[0]);
    }
  }
}

void Ahadic_Reweighting::KTSelectionReweighting(
    const bool accepted, const std::vector<double>& probs) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from KT_Selector for each accept/reject trial of the kt
  // selection. NOTE: these weights are computed but never committed --
  // see the TODO in KT_Selector::operator() -- so KT_0 variations are
  // incomplete for the time being.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    if (accepted) {
      m_tmp_kt_weights[ivar] *= probs[ivar] / probs[0];
    } else {
      m_tmp_kt_weights[ivar] *= (1. - probs[ivar]) / (1. - probs[0]);
    }
  }
}

void Ahadic_Reweighting::SoftClusterReweighting(
    const std::vector<double>& weights, const std::vector<double>& totweights) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Soft_Cluster_Handler::DecayWeight for the selected hadron
  // pair. The selection probability is p = weight / totweight with the
  // variant double-transition tables; weight *= p_var / p_nom. These
  // weights are committed directly into the event weights.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  const double prob_nom = weights[0] / totweights[0];
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    // TODO: figure out why this is zero from time to time
    const double ratio = (weights[ivar] / totweights[ivar]) / prob_nom;
    if (!std::isnan(ratio)) {
      m_soft_weights[ivar] *= ratio;
    }
  }
}

void Ahadic_Reweighting::ApplyVariationWeights(ATOOLS::Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Compute and apply the variation weights of one hadronization call.
  // The total weight is:
  // w_total = w_cluster * w_gluon * w_flavour * w_kt * w_soft
  // capped at m_max_reweight_factor, and multiplied into the soft-physics
  // variations of the event's Weights_Map. Hadronization may run more than
  // once per event; each call multiplies its weights into the same variations, 
  // so the event weight is the product over all calls.
  ///////////////////////////////////////////////////////////////////////////
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    double w_total = 1.;
    w_total *= m_cluster_weights[ivar];
    w_total *= m_gluon_weights[ivar];
    w_total *= m_flavour_weights[ivar];
    w_total *= m_kt_weights[ivar];
    w_total *= m_soft_weights[ivar];
    if (std::isnan(w_total)) {
      msg_Error() << METHOD << ": NaN variation weight, resetting to 1\n";
      w_total = 1.0;
    }
    if (m_max_reweight_factor > 0. && w_total > m_max_reweight_factor) {
      w_total = m_max_reweight_factor;
    }
    m_variation_weights[ivar] = w_total;
  }
  if (blob != NULL) {
    auto wgtmap = (*blob)["WeightsMap"]->Get<Weights_Map>();
    CombineSoftPhysicsVariations(wgtmap, m_variation_weights);
    blob->AddData("WeightsMap", new Blob_Data<Weights_Map>(wgtmap));
  }
  ResetEvent();
}
