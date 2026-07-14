#include "AHADIC++/Tools/Hadronisation_Reweighting.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <algorithm>
#include <cmath>
#include "ATOOLS/Org/MyStrStream.H" // OUTPUT
#include "ATOOLS/Math/Histogram.H"  // OUTPUT
#include "ATOOLS/Org/Shell_Tools.H" // OUTPUT
#include <sstream>                  // OUTPUT

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Hadronisation_Reweighting::Hadronisation_Reweighting() :
  m_n_variations(1),
  m_max_reweight_factor(-1.),
  m_reweight_max_nsplit(-1),

  // OUTPUT
  m_output_mode(0),
  m_reweighting_output(false),
  m_event_weights_applied(false),
  m_total_events(0)
{}

Hadronisation_Reweighting::~Hadronisation_Reweighting() {
  // OUTPUT
  PrintVariationStatistics();
  if (m_hadronisation_weight_file.is_open()) {
    m_hadronisation_weight_file.close();
  }
  if (m_output_mode == 2) WriteHistograms();
}

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

  // OUTPUT
  ResetStats();
  auto s = Settings::GetMainSettings()["AHADIC"];
  const int reweighting_output = s["REWEIGHTING_OUTPUT"].SetDefault(0).Get<int>();
  m_output_mode        = reweighting_output;
  m_reweighting_output = reweighting_output > 0;
  m_cutoff_count.resize(m_n_variations, 0);
  m_sum_weights.resize(m_n_variations, 0.0);
  m_sum_weights_squared.resize(m_n_variations, 0.0);
  m_total_events = 0;
  if (m_output_mode == 1) {
    m_hadronisation_weight_file.open("hadronisation_weights.dat");
    if (m_hadronisation_weight_file.is_open()) {
      m_hadronisation_weight_file
        << "# W <event> <ivar> <w_flavour> <w_gluon> <w_cluster> <w_soft> <w_kt> <w_prompt> <w_total>\n"
        << "# G <event> <z>\n"
        << "# C <event> <nsplit> <type1> <z1> <type2> <z2> <M>\n"
        << "# K <event> <site:0=gluon,1=cluster,2=soft> <kt> <ktmax>\n"
        << "# F <event> <kf> <Emax>\n"
        << "# S <event> <M> <kf1> <kf2>\n"
        << "# T <event> <M> <kf>\n"
        << "# P <event> <M> <kf1> <kf2>\n"
        << "# E <event> <n_primary> <n_gluon> <n_fission> <n_soft>\n";
      m_hadronisation_weight_file << std::scientific << std::setprecision(10);
    }
  }
  else if (m_output_mode == 2) {
    BookHistograms();
  }
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
  // OUTPUT
  m_tmp_flavour_records.clear();
  m_tmp_gluon_records.clear();
  m_tmp_cluster_records.clear();
  m_call_records.clear();
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
      m_cutoff_count[ivar]++; // OUTPUT
    }
    m_variation_weights[ivar] = w_total;
  }
  if (blob != NULL) {
    auto wgtmap = (*blob)["WeightsMap"]->Get<Weights_Map>();
    CombineSoftPhysicsVariations(wgtmap, m_variation_weights);
    blob->AddData("WeightsMap", new Blob_Data<Weights_Map>(wgtmap));
    AccumulateEventStatistics(); // OUTPUT
  }
  ResetEvent();
}

///////////////////////////// OUTPUT AND STATISTICS /////////////////////////////

void Hadronisation_Reweighting::ResetStats() {
  m_event_weights_applied = false;
  m_event_variation_weights.assign(m_n_variations, 1.);
  m_event_flavour_weights.assign(m_n_variations, 1.);
  m_event_gluon_weights.assign(m_n_variations, 1.);
  m_event_cluster_weights.assign(m_n_variations, 1.);
  m_event_soft_weights.assign(m_n_variations, 1.);
  m_event_kt_weights.assign(m_n_variations, 1.);
  m_event_promptdecay_weights.assign(m_n_variations, 1.);
  m_event_records.clear();
}

void Hadronisation_Reweighting::AccumulateEventStatistics() {
  ///////////////////////////////////////////////////////////////////////////
  // Multiply the weights of one hadronization call into the per-event
  // weights; hadronization may run more than once per event.
  ///////////////////////////////////////////////////////////////////////////
  m_event_weights_applied = true;
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    m_event_flavour_weights[ivar]     *= m_flavour_weights[ivar];
    m_event_gluon_weights[ivar]       *= m_gluon_weights[ivar];
    m_event_cluster_weights[ivar]     *= m_cluster_weights[ivar];
    m_event_soft_weights[ivar]        *= m_soft_weights[ivar];
    m_event_kt_weights[ivar]          *= m_kt_weights[ivar];
    m_event_promptdecay_weights[ivar] *= m_promptdecay_weights[ivar];
    m_event_variation_weights[ivar]   *= m_variation_weights[ivar];
  }
  m_event_records.insert(m_event_records.end(),
                         m_call_records.begin(), m_call_records.end());
  m_call_records.clear();
}

void Hadronisation_Reweighting::WriteEventStatistics() {
  ///////////////////////////////////////////////////////////////////////////
  // Flush the per-event weights into the statistics; called once per event
  // when the event phases are cleaned up.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_event_weights_applied) return;
  m_total_events++;
  if (m_hadronisation_weight_file.is_open()) {
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      m_hadronisation_weight_file << "W " << m_total_events << " " << ivar
                    << " " << m_event_flavour_weights[ivar]
                    << " " << m_event_gluon_weights[ivar]
                    << " " << m_event_cluster_weights[ivar]
                    << " " << m_event_soft_weights[ivar]
                    << " " << m_event_kt_weights[ivar]
                    << " " << m_event_promptdecay_weights[ivar]
                    << " " << m_event_variation_weights[ivar] << "\n";
    }
    size_t n_primary(0), n_gluon(0), n_fission(0), n_soft(0);
    for (const auto& record : m_event_records) {
      switch (record.first) {
      case 'P': ++n_primary; break;
      case 'G': ++n_gluon;   break;
      case 'C': ++n_fission; break;
      case 'S': ++n_soft;    break;
      default: break;
      }
      m_hadronisation_weight_file << record.first << " " << m_total_events
                    << record.second << "\n";
    }
    m_hadronisation_weight_file << "E " << m_total_events << " " << n_primary << " "
                  << n_gluon << " " << n_fission << " " << n_soft << "\n";
  }
  if (m_output_mode == 2) FillHistograms();
  for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
    const double w = m_event_variation_weights[ivar];
    m_sum_weights[ivar] += w;
    m_sum_weights_squared[ivar] += w * w;
  }
  ResetStats();
}

namespace {
  const long int s_flavour_codes[] = {
    3303, 3201, 3101, 3203, 3103, 2101, 2103, 2203, 1103, 3, 2, 1
  };
  const int s_n_flavour = int(sizeof(s_flavour_codes)/sizeof(long int));
  const long int s_hadron_codes[] = {
    22, 111, 211, 221, 331, 113, 213, 223, 333, 311, 321, 313, 323,
    215, 225, 315, 325, 411, 421, 431, 413, 423, 433, 511, 521, 531,
    513, 523, 533, 2212, 2112, 1114, 2114, 2214, 2224, 3122, 3112,
    3212, 3222, 3114, 3214, 3224, 3312, 3322, 3314, 3324, 3334, 4122, 5122
  };
  const int s_n_hadron = int(sizeof(s_hadron_codes)/sizeof(long int));

  inline long int AbsKF(const long int kf) { return kf < 0 ? -kf : kf; }
  int FlavourCategory(const long int kf) {
    const long int akf = AbsKF(kf);
    for (int i=0; i<s_n_flavour; ++i)
      if (s_flavour_codes[i]==akf) return i;
    return -1;
  }
  int HadronCategory(const long int kf) {
    const long int akf = AbsKF(kf);
    for (int i=0; i<s_n_hadron; ++i)
      if (s_hadron_codes[i]==akf) return i;
    return s_n_hadron;
  }
}

void Hadronisation_Reweighting::BookHistograms() {
  ///////////////////////////////////////////////////////////////////////////
  // Book one nominal and (per variation) one reweighted histogram for every
  // observable, with fixed, hardcoded binning (type 1 = linear + error).
  ///////////////////////////////////////////////////////////////////////////
  m_hist_path = "AHADIC_Histograms";
  MakeDir(m_hist_path);
  auto book = [this](const std::string& key,
                     const double lo, const double hi, const int nbins) {
    m_hist_nominal[key] = new Histogram(1, lo, hi, nbins, key);
    std::vector<Histogram*> vec(m_n_variations, (Histogram*)NULL);
    for (size_t ivar=1; ivar<m_n_variations; ++ivar)
      vec[ivar] = new Histogram(1, lo, hi, nbins, key);
    m_hist_reweighted[key] = vec;
  };
  book("gluon_z",           0.,  1., 80);
  book("cluster_z_light",   0.,  1., 80);
  book("cluster_z_leading", 0.,  1., 80);
  book("cluster_z_diquark", 0.,  1., 80);
  book("cluster_z_beam",    0.,  1., 80);
  book("kt_gluon",          0.,  3., 80);
  book("kt_cluster",        0.,  3., 80);
  book("kt_soft",           0.,  3., 80);
  book("flavour_pop",     -0.5, s_n_flavour-0.5, s_n_flavour);
  book("soft_hadron",     -0.5, s_n_hadron+0.5,  s_n_hadron+1);
  book("transition_hadron", -0.5, s_n_hadron+0.5, s_n_hadron+1);
  book("soft_cluster_mass",       0.,  5., 80);
  book("transition_cluster_mass", 0.,  5., 80);
  book("cluster_mass",            0.,  5., 80);
  book("primary_cluster_mass",    0., 20., 80);
  book("cluster_nsplit",       -0.5, 30.5, 31);
  book("n_primary_clusters",   -0.5, 100.5, 101);
  book("n_gluon_splittings",   -0.5, 100.5, 101);
  book("n_cluster_splittings", -0.5, 100.5, 101);
  book("n_soft_decays",        -0.5, 100.5, 101);
}

void Hadronisation_Reweighting::FillObs(const std::string& key,
                                        const double value) {
  ///////////////////////////////////////////////////////////////////////////
  // Fill one observable value into the nominal histogram (weight 1) and into
  // each reweighted histogram with the finalised per-event variation weight.
  ///////////////////////////////////////////////////////////////////////////
  std::map<std::string, Histogram*>::iterator nit = m_hist_nominal.find(key);
  if (nit == m_hist_nominal.end()) return;
  nit->second->Insert(value, 1.0);
  std::vector<Histogram*>& vec = m_hist_reweighted[key];
  for (size_t ivar=1; ivar<m_n_variations; ++ivar)
    vec[ivar]->Insert(value, m_event_variation_weights[ivar]);
}

void Hadronisation_Reweighting::FillHistograms() {
  ///////////////////////////////////////////////////////////////////////////
  // Replay the per-event records into the histograms; the reweighting weight
  // is the per-event total, applied uniformly to every record of the event.
  ///////////////////////////////////////////////////////////////////////////
  static const char* cluster_class[4] = {"light", "leading", "diquark", "beam"};
  size_t n_primary(0), n_gluon(0), n_fission(0), n_soft(0);
  for (const auto& record : m_event_records) {
    std::istringstream str(record.second);
    switch (record.first) {
    case 'G': {
      double z; str >> z;
      FillObs("gluon_z", z);
      ++n_gluon;
      break;
    }
    case 'C': {
      int nsplit, t1, t2; double z1, z2, mass;
      str >> nsplit >> t1 >> z1 >> t2 >> z2 >> mass;
      FillObs("cluster_nsplit", double(nsplit));
      FillObs("cluster_mass", mass);
      if (t1>=0 && t1<4)
        FillObs(std::string("cluster_z_")+cluster_class[t1], z1);
      if (t2>=0 && t2<4)
        FillObs(std::string("cluster_z_")+cluster_class[t2], z2);
      ++n_fission;
      break;
    }
    case 'K': {
      int site; double kt, ktmax;
      str >> site >> kt >> ktmax;
      if      (site==0) FillObs("kt_gluon", kt);
      else if (site==1) FillObs("kt_cluster", kt);
      else if (site==2) FillObs("kt_soft", kt);
      break;
    }
    case 'F': {
      long int kf; double Emax;
      str >> kf >> Emax;
      const int cat = FlavourCategory(kf);
      if (cat>=0) FillObs("flavour_pop", double(cat));
      break;
    }
    case 'S': {
      double mass; long int kf1, kf2;
      str >> mass >> kf1 >> kf2;
      FillObs("soft_cluster_mass", mass);
      FillObs("soft_hadron", double(HadronCategory(kf1)));
      FillObs("soft_hadron", double(HadronCategory(kf2)));
      ++n_soft;
      break;
    }
    case 'T': {
      double mass; long int kf;
      str >> mass >> kf;
      FillObs("transition_cluster_mass", mass);
      FillObs("transition_hadron", double(HadronCategory(kf)));
      break;
    }
    case 'P': {
      double mass; long int kf1, kf2;
      str >> mass >> kf1 >> kf2;
      FillObs("primary_cluster_mass", mass);
      ++n_primary;
      break;
    }
    default: break;
    }
  }
  FillObs("n_primary_clusters",   double(n_primary));
  FillObs("n_gluon_splittings",   double(n_gluon));
  FillObs("n_cluster_splittings", double(n_fission));
  FillObs("n_soft_decays",        double(n_soft));
}

void Hadronisation_Reweighting::WriteHistograms() {
  ///////////////////////////////////////////////////////////////////////////
  // Write out the raw (un-normalised) histograms, one file per histogram, so
  // that sub-runs can be merged by bin-wise addition and normalised at plot
  // time.
  ///////////////////////////////////////////////////////////////////////////
  for (auto& it : m_hist_nominal) {
    Histogram* h = it.second;
    h->Output(m_hist_path + "/" + it.first + "_nominal.dat");
    delete h;
  }
  m_hist_nominal.clear();
  for (auto& it : m_hist_reweighted) {
    std::vector<Histogram*>& vec = it.second;
    for (size_t ivar=1; ivar<vec.size(); ++ivar) {
      if (vec[ivar]==NULL) continue;
      vec[ivar]->Output(m_hist_path + "/" + it.first +
                        "_reweighted_v" + std::to_string(ivar) + ".dat");
      delete vec[ivar];
    }
  }
  m_hist_reweighted.clear();
}

void Hadronisation_Reweighting::CommitTmpRecords(Record_Vector& tmp_records) {
  ///////////////////////////////////////////////////////////////////////////
  // Commit the records of a successful splitting into the call records.
  ///////////////////////////////////////////////////////////////////////////
  if (tmp_records.empty()) return;
  m_call_records.insert(m_call_records.end(),
                        tmp_records.begin(), tmp_records.end());
  tmp_records.clear();
}

void Hadronisation_Reweighting::DiscardTmpRecords(Record_Vector& tmp_records) {
  ///////////////////////////////////////////////////////////////////////////
  // Discard the records of a failed splitting attempt.
  ///////////////////////////////////////////////////////////////////////////
  tmp_records.clear();
}

void Hadronisation_Reweighting::RecordFlavourPop(const Flavour& flav,
                                          const double& Emax) {
  ///////////////////////////////////////////////////////////////////////////
  // F record from Flavour_Selector: the popped flavour and the Emax cut of
  // the event-dependent normalisation. Committed with the flavour weights,
  // i.e. including pops of failed attempts of a successful splitting.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_reweighting_output) return;
  MyStrStream str;
  str << " " << (long int)(flav)
      << " " << std::scientific << std::setprecision(8) << Emax;
  m_tmp_flavour_records.emplace_back('F', str.str());
}

void Hadronisation_Reweighting::RecordGluonZ(const double& z) {
  ///////////////////////////////////////////////////////////////////////////
  // G record from Gluon_Splitter: the realized z of a gluon splitting.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_reweighting_output) return;
  MyStrStream str;
  str << " " << std::scientific << std::setprecision(8) << z;
  m_tmp_gluon_records.emplace_back('G', str.str());
}

void Hadronisation_Reweighting::RecordClusterZ(const int nsplit,
                                        const int type1, const double& z1,
                                        const int type2, const double& z2,
                                        const double& mass) {
  ///////////////////////////////////////////////////////////////////////////
  // C record from Cluster_Splitter: chain index and the realized z per
  // cluster side with its flavour class (0=light, 1=leading, 2=diquark,
  // 3=beam), plus the mass of the decaying cluster.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_reweighting_output) return;
  MyStrStream str;
  str << " " << nsplit << " " << std::scientific << std::setprecision(8)
      << type1 << " " << z1 << " " << type2 << " " << z2 << " " << mass;
  m_tmp_cluster_records.emplace_back('C', str.str());
}

void Hadronisation_Reweighting::RecordKT(const int site, const double& kt,
                                  const double& ktmax) {
  ///////////////////////////////////////////////////////////////////////////
  // K record: a kt draw with its event-dependent cut. Splitter draws
  // (site 0/1) share the fate of the splitting attempt; soft-cluster draws
  // (site 2) are committed directly. Untagged sites (rescue paths) and
  // failed draws are not recorded.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_reweighting_output) return;
  if (site < 0 || kt < 0.) return;
  MyStrStream str;
  str << " " << site << " " << std::scientific << std::setprecision(8)
      << kt << " " << ktmax;
  if      (site == 0) m_tmp_gluon_records.emplace_back('K', str.str());
  else if (site == 1) m_tmp_cluster_records.emplace_back('K', str.str());
  else                m_call_records.emplace_back('K', str.str());
}

void Hadronisation_Reweighting::RecordSoftDecay(const double& mass,
                                         const Flavour& had1,
                                         const Flavour& had2) {
  ///////////////////////////////////////////////////////////////////////////
  // S record from Soft_Cluster_Handler::DecayWeight: cluster mass and the
  // selected hadron pair. Committed directly, mirroring the soft weights
  // (also for decay attempts that fail later).
  ///////////////////////////////////////////////////////////////////////////
  if (!m_reweighting_output) return;
  MyStrStream str;
  str << " " << std::scientific << std::setprecision(8) << mass
      << " " << (long int)(had1) << " " << (long int)(had2);
  m_call_records.emplace_back('S', str.str());
}

void Hadronisation_Reweighting::RecordTransition(const double& mass,
                                          const Flavour& had) {
  ///////////////////////////////////////////////////////////////////////////
  // T record from Soft_Cluster_Handler::RadiationWeight: cluster mass and
  // the selected single-hadron transition. Committed directly, mirroring
  // the soft weights (also for transitions that are not realized later).
  ///////////////////////////////////////////////////////////////////////////
  if (!m_reweighting_output) return;
  MyStrStream str;
  str << " " << std::scientific << std::setprecision(8) << mass
      << " " << (long int)(had);
  m_call_records.emplace_back('T', str.str());
}

void Hadronisation_Reweighting::RecordPrimaryCluster(const double& mass,
                                              const Flavour& fl1,
                                              const Flavour& fl2,
                                              const bool direct) {
  ///////////////////////////////////////////////////////////////////////////
  // P record: mass and constituent flavours of a primary cluster. From
  // Gluon_Splitter the record shares the fate of the splitting attempt;
  // directly formed clusters (Gluon_Decayer::Trivial) commit immediately.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_reweighting_output) return;
  MyStrStream str;
  str << " " << std::scientific << std::setprecision(8) << mass
      << " " << (long int)(fl1) << " " << (long int)(fl2);
  if (direct) m_call_records.emplace_back('P', str.str());
  else        m_tmp_gluon_records.emplace_back('P', str.str());
}

void Hadronisation_Reweighting::PrintVariationStatistics() {
  if (m_n_variations <= 1 || m_total_events == 0) return;

  const std::string title = "AHADIC Reweighting Statistics (events: " +
                            ToString<size_t>(m_total_events) + ")";

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
