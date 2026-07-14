#include "AHADIC++/Tools/Hadronisation_Reweighting.H"
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

Hadronisation_Reweighting::Hadronisation_Reweighting() :
  m_n_variations(1),
  m_max_reweight_factor(-1.),
  m_reweight_max_nsplit(-1)
{}

Hadronisation_Reweighting::~Hadronisation_Reweighting() {}

void Hadronisation_Reweighting::Initialize() {
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
    "kT_max",
    "mass_exponent",
    "prompt_decay_exponent",
    "Singlet_Suppression",
    "Mixing_Angle_0+",
    "Mixing_Angle_1-",
    "Mixing_Angle_2+",
    "Multiplet_Meson_R0L0S0",
    "Multiplet_Meson_R0L0S1",
    "Multiplet_Meson_R0L0S2",
    "Multiplet_Meson_R0L1S0",
    "Multiplet_Meson_R0L1S1",
    "Multiplet_Meson_R0L2S2",
    "Multiplet_Baryon_R0L0S1/2",
    "Multiplet_Baryon_R1L0S1/2",
    "Multiplet_Baryon_R2L0S1/2",
    "Multiplet_Baryon_R1_1L0S1/2",
    "Multiplet_Baryon_R0L0S3/2",
    "eta_modifier",
    "eta_prime_modifier",
    "Singlet_Baryon_modifier",
    "CharmBaryon_Enhancement",
    "BeautyBaryon_Enhancement",
    "CharmStrange_Enhancement",
    "BeautyStrange_Enhancement",
    "BeautyCharm_Enhancement",
    "Strange_fraction",
    "Baryon_fraction",
    "P_qs_by_P_qq",
    "P_ss_by_P_qq",
    "P_di_1_by_P_di_0"
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
  if (m_n_variations > 1) CheckVariationGuards();
  for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
    const double strange = m_variation_vectors["Strange_fraction"][ivar];
    m_variation_vectors["P_qs_by_P_qq"][ivar] *= strange;
    m_variation_vectors["P_ss_by_P_qq"][ivar] *= sqr(strange);
  }

  m_max_reweight_factor = hadpars->Get("max_reweight_factor");
  m_reweight_max_nsplit = hadpars->Switch("reweight_max_nsplit");
  ResetEvent();
}

void Hadronisation_Reweighting::CheckVariationGuards() {
  ///////////////////////////////////////////////////////////////////////////
  // Guards on the parameter variations, checked once at initialisation on
  // the raw parameter vectors. The composite table weights are guarded by
  // CheckTransitionGuard / CheckPoppingGuard / CheckComponentGuard below,
  // called from Single_Transitions::FillMap, Double_Transitions::FillMap
  // and Multiplet_Constructor::ComponentActive during table construction.
  ///////////////////////////////////////////////////////////////////////////
  const std::vector<double>& ktmax = m_variation_vectors["kT_max"];
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    if (ktmax[ivar] > ktmax[0]) {
      THROW(fatal_error, std::string("Reweighting of AHADIC only possible for downward variations of PT_MAX.\n")
                        + "Found PT_MAX variation " + std::to_string(ivar) + " = "
                        + std::to_string(ktmax[ivar]) + " > "
                        + std::to_string(ktmax[0]) + " = PT_MAX nominal.\n"
                        + "Please adjust your parameter variation settings.");
    }
  }
  for (const auto& kv : m_variation_vectors) {
    if (kv.first.rfind("Multiplet_", 0) != 0) continue;
    const std::vector<double>& wts = kv.second;
    if (wts[0] < 1.e-6) {
      for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
        if (wts[ivar] >= 1.e-6) {
          THROW(fatal_error, std::string("Reweighting of AHADIC not possible for multiplets switched off in the nominal run.\n")
                            + "Found " + kv.first + " variation " + std::to_string(ivar) + " = "
                            + std::to_string(wts[ivar]) + " with nominal weight "
                            + std::to_string(wts[0]) + ".\n"
                            + "Please adjust your parameter variation settings.");
        }
      }
    }
    else {
      for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
        if (wts[ivar] < 1.e-6) {
          THROW(fatal_error, std::string("Reweighting of AHADIC not possible for multiplets switched off by a variation.\n")
                            + "Found " + kv.first + " variation " + std::to_string(ivar) + " = "
                            + std::to_string(wts[ivar]) + " with nominal weight "
                            + std::to_string(wts[0]) + ".\n"
                            + "Switching a multiplet off changes the transition tables and their mass thresholds;\n"
                            + "use a small non-zero weight instead.\n"
                            + "Please adjust your parameter variation settings.");
        }
      }
    }
  }
  static const std::vector<std::pair<std::string, std::string>> popping_keys = {
    {"Strange_fraction", "STRANGE_FRACTION"},
    {"Baryon_fraction",  "BARYON_FRACTION"},
    {"P_qs_by_P_qq",     "P_QS_by_P_QQ_norm"},
    {"P_ss_by_P_qq",     "P_SS_by_P_QQ_norm"},
    {"P_di_1_by_P_di_0", "P_QQ1_by_P_QQ0"},
  };
  for (const auto& kv : popping_keys) {
    const std::vector<double>& vals = m_variation_vectors[kv.first];
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      if (vals[0] <= 0. && vals[ivar] > 0.) {
        THROW(fatal_error, std::string("Reweighting of AHADIC not possible for flavour popping switched off in the nominal run.\n")
                          + "Found " + kv.second + " variation " + std::to_string(ivar) + " = "
                          + std::to_string(vals[ivar]) + " with nominal value "
                          + std::to_string(vals[0]) + ".\n"
                          + "Please adjust your parameter variation settings.");
      }
    }
  }
  const std::vector<double>& zeta = m_variation_vectors["prompt_decay_exponent"];
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    if ((zeta[0] > 0.) != (zeta[ivar] > 0.)) {
      THROW(fatal_error, std::string("Reweighting of AHADIC prompt-decay exponent only possible within the stochastic mode (zeta > 0).\n")
                        + "Found PROMPT_DECAY_EXPONENT variation " + std::to_string(ivar) + " = "
                        + std::to_string(zeta[ivar]) + " with nominal "
                        + std::to_string(zeta[0]) + ".\n"
                        + "Please adjust your parameter variation settings.");
    }
  }
}

void Hadronisation_Reweighting::CheckTransitionGuard(
    const ATOOLS::Flavour& hadron, const std::vector<double>& weights,
    const double& threshold) const {
  ///////////////////////////////////////////////////////////////////////////
  // Guard called from Single_Transitions::FillMap: a hadron's table
  // membership must not differ between the nominal run and any variation;
  // a dedicated run at the varied value would prune or add the hadron, shifting
  // the transition and decay thresholds derived from the tables' mass extremes.
  ///////////////////////////////////////////////////////////////////////////
  if (weights[0] < threshold) {
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      if (weights[ivar] >= threshold) {
        THROW(fatal_error, std::string("Reweighting of AHADIC not possible for hadron transitions ")
                    + "absent from the nominal transition tables.\n"
                    + "Variation " + std::to_string(ivar) + " is activating the hadron "
                    + ToString(hadron) + ", deactivated in the nominal run.\n"
                    + "Please adjust your parameter variation settings.");
      }
    }
    return;
  }
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    if (weights[ivar] < threshold) {
      THROW(fatal_error, std::string("Reweighting of AHADIC not possible for hadron transitions ")
                  + "switched off by a variation.\n"
                  + "Variation " + std::to_string(ivar) + " is deactivating the hadron "
                  + ToString(hadron) + ", present in the nominal transition tables.\n"
                  + "Use a small non-zero weight instead.\n"
                  + "Please adjust your parameter variation settings.");
    }
  }
}

void Hadronisation_Reweighting::CheckPoppingGuard(
    const ATOOLS::Flavour& popped, const std::vector<double>& weights,
    const double& threshold) const {
  ///////////////////////////////////////////////////////////////////////////
  // Guard called from Double_Transitions::FillMap: a popped flavour's
  // constituent weight must stay on the same side of the table threshold in
  // all variations; a dedicated run at the varied value would prune or
  // add all decay channels with this popped flavour, shifting the decay
  // thresholds.
  ///////////////////////////////////////////////////////////////////////////
  if (weights[0] < threshold) {
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      if (weights[ivar] >= threshold) {
        THROW(fatal_error, std::string("Reweighting of AHADIC not possible for flavour popping ")
              + "absent from the nominal run.\n"
              + "Variation " + std::to_string(ivar) + " is activating popping of "
              + ToString(popped) + ", deactivated in the nominal run.\n"
              + "Please adjust your parameter variation settings.");
      }
    }
    return;
  }
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    if (weights[ivar] < threshold) {
      THROW(fatal_error, std::string("Reweighting of AHADIC not possible for flavour popping ")
            + "switched off by a variation.\n"
            + "Variation " + std::to_string(ivar) + " is deactivating popping of "
            + ToString(popped) + ", active in the nominal run.\n"
            + "Use a small non-zero value instead.\n"
            + "Please adjust your parameter variation settings.");
    }
  }
}

void Hadronisation_Reweighting::CheckComponentGuard(
    const ATOOLS::Flavour& hadron, const std::vector<double>& amps,
    const double& threshold) const {
  ///////////////////////////////////////////////////////////////////////////
  // Guard called from Multiplet_Constructor::ComponentActive: a wave
  // function component active in the nominal run must stay active in all
  // variations and vice versa; a dedicated run at the varied mixing angle
  // would drop or add the component, pruning the hadron from the
  // corresponding transition lists.
  ///////////////////////////////////////////////////////////////////////////
  if (dabs(amps[0]) > threshold) {
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      if (dabs(amps[ivar]) <= threshold) {
        THROW(fatal_error, std::string("Reweighting of AHADIC not possible for wave function components ")
                          + "switched off by a variation.\n"
                          + "The mixing-angle variation " + std::to_string(ivar) + " is deactivating a component of "
                          + ToString(hadron) + ", active in the nominal run.\n"
                          + "Please adjust your parameter variation settings.");
      }
    }
    return;
  }
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    if (dabs(amps[ivar]) > threshold) {
      THROW(fatal_error, std::string("Reweighting of AHADIC not possible for wave function components ")
                        + "absent from the nominal run.\n"
                        + "The mixing-angle variation " + std::to_string(ivar) + " is activating a component of "
                        + ToString(hadron) + ", deactivated in the nominal run.\n"
                        + "Please adjust your parameter variation settings.");
    }
  }
}

std::vector<double>
Hadronisation_Reweighting::GetVariationVector(const std::string& keyword) const {
  ///////////////////////////////////////////////////////////////////////////
  // Access to the padded parameter-variation vectors for the tools that
  // compute the per-variation selection probabilities.
  ///////////////////////////////////////////////////////////////////////////
  auto piter = m_variation_vectors.find(keyword);
  if (piter != m_variation_vectors.end()) return piter->second;
  THROW(fatal_error, "Keyword " + keyword +
                     " not found in Hadronisation_Reweighting vector map.");
}

void Hadronisation_Reweighting::ResetEvent() {
  ///////////////////////////////////////////////////////////////////////////
  // Reset all variation weights to 1.0 for a new event.
  ///////////////////////////////////////////////////////////////////////////
  m_variation_weights.assign(m_n_variations, 1.);
  m_flavour_weights.assign(m_n_variations, 1.);
  m_gluon_weights.assign(m_n_variations, 1.);
  m_cluster_weights.assign(m_n_variations, 1.);
  m_soft_weights.assign(m_n_variations, 1.);
  m_kt_weights.assign(m_n_variations, 1.);
  m_promptdecay_weights.assign(m_n_variations, 1.);
  m_tmp_flavour_weights.assign(m_n_variations, 1.);
  m_tmp_gluon_weights.assign(m_n_variations, 1.);
  m_tmp_cluster_weights.assign(m_n_variations, 1.);
}

void Hadronisation_Reweighting::AcceptTmpWeights(std::vector<double>& weights,
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

void Hadronisation_Reweighting::ResetTmpWeights(std::vector<double>& tmp_weights) {
  ///////////////////////////////////////////////////////////////////////////
  // Discard the weights of a failed splitting attempt.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  std::fill(tmp_weights.begin(), tmp_weights.end(), 1.);
}

void Hadronisation_Reweighting::FlavourSelectionReweighting(
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

void Hadronisation_Reweighting::GluonSplittingReweighting(
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

void Hadronisation_Reweighting::ClusterSplittingReweighting(
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
      if (std::isfinite(ratio)) {
        m_tmp_cluster_weights[ivar] *= ratio;
      }
    } else {
      m_tmp_cluster_weights[ivar] *= (1. - probs[ivar]) / (1. - probs[0]);
    }
  }
}

void Hadronisation_Reweighting::KTSelectionReweighting(
    const std::vector<double>& probs) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from KT_Selector for each completed kt draw. The probs are the
  // normalised truncated-Gaussian densities of the realized kt, with a
  // variation's prob set to zero if kt lies above its tightened PT_MAX cut,
  // so that the accept/reject trials integrate out and rejected kt play no
  // role; weight *= p_var / p_nom. Committed directly into the event
  // weights, also when the enclosing splitting attempt fails later.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    m_kt_weights[ivar] *= probs[ivar] / probs[0];
  }
}

void Hadronisation_Reweighting::SoftClusterReweighting(
    const std::vector<double>& weights, const std::vector<double>& totweights) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Soft_Cluster_Handler for the selected hadron (pair). The
  // selection probability is p = weight / totweight with the variant
  // transition tables; weight *= p_var / p_nom. These weights are
  // committed directly into the event weights.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  const double prob_nom = weights[0] / totweights[0];
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    const double ratio = (weights[ivar] / totweights[ivar]) / prob_nom;
    if (std::isfinite(ratio)) {
      m_soft_weights[ivar] *= ratio;
    }
  }
}

void Hadronisation_Reweighting::PromptDecayReweighting(
    const bool decayed, const std::vector<double>& probs) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Soft_Cluster_Handler::MustPromptDecay for the Bernoulli
  // prompt-decay decision of a cluster, only taken in the stochastic mode
  // (zeta > 0); probs are the per-variation decay probabilities.
  // For a decaying cluster:   weight *= p_var / p_nom
  // For a splitting cluster:  weight *= (1 - p_var) / (1 - p_nom)
  // Committed directly into the event weights.
  ///////////////////////////////////////////////////////////////////////////
  if (m_n_variations <= 1) return;
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    m_promptdecay_weights[ivar] *= decayed ?
      probs[ivar] / probs[0] : (1. - probs[ivar]) / (1. - probs[0]);
  }
}

void Hadronisation_Reweighting::ApplyVariationWeights(ATOOLS::Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Compute and apply the variation weights of one hadronization call.
  // The total weight is:
  // w_total = w_cluster * w_gluon * w_flavour * w_kt * w_soft * w_prompt
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
    w_total *= m_promptdecay_weights[ivar];
    if (!std::isfinite(w_total)) {
      msg_Error() << METHOD << ": non-finite variation weight, resetting to 1\n";
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
