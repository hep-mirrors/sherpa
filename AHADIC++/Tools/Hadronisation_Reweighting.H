#ifndef AHADIC_Tools_Hadronisation_Reweighting_H
#define AHADIC_Tools_Hadronisation_Reweighting_H

#include "ATOOLS/Phys/Blob.H"
#include <vector>
#include <map>
#include <string>

namespace AHADIC {

  class Hadronisation_Reweighting {
  private:
    size_t m_n_variations;
    double m_max_reweight_factor;
    int    m_reweight_max_nsplit;
    std::map<std::string, std::vector<double>> m_variation_vectors;

    std::vector<double> m_variation_weights;
    // Per-channel event weights, accumulated over all splittings of an event:
    // flavour popping (Flavour_Selector), gluon splitting z (Gluon_Splitter),
    // cluster splitting z (Cluster_Splitter), soft cluster decays and
    // prompt-decay decisions (Soft_Cluster_Handler), and transverse
    // momentum (KT_Selector).
    std::vector<double> m_flavour_weights;
    std::vector<double> m_gluon_weights;
    std::vector<double> m_cluster_weights;
    std::vector<double> m_soft_weights;
    std::vector<double> m_kt_weights;
    std::vector<double> m_promptdecay_weights;
    // Per-splitting weights, committed into the event weights when a
    // splitting succeeds and discarded when an attempt fails.
    std::vector<double> m_tmp_flavour_weights;
    std::vector<double> m_tmp_gluon_weights;
    std::vector<double> m_tmp_cluster_weights;

    void CheckVariationGuards();
    void AcceptTmpWeights(std::vector<double>& weights,
                       std::vector<double>& tmp_weights);
    void ResetTmpWeights(std::vector<double>& tmp_weights);
  public:
    Hadronisation_Reweighting();
    ~Hadronisation_Reweighting();

    void Initialize();
    void ResetEvent();

    inline bool   Active() const { return m_n_variations > 1; }
    inline size_t NumberOfVariations() const { return m_n_variations; }
    std::vector<double> GetVariationVector(const std::string& keyword) const;

    void CheckTransitionGuard(const ATOOLS::Flavour& hadron,
                              const std::vector<double>& weights,
                              const double& threshold) const;
    void CheckPoppingGuard(const ATOOLS::Flavour& popped,
                           const std::vector<double>& weights,
                           const double& threshold) const;
    void CheckComponentGuard(const ATOOLS::Flavour& hadron,
                             const std::vector<double>& amps,
                             const double& threshold) const;

    // Flavour popping in the shared Flavour_Selector
    void FlavourSelectionReweighting(const std::vector<double>& popweights,
                                     const std::vector<double>& norms);
    inline void AcceptFlavourSelectionWeights() {
      AcceptTmpWeights(m_flavour_weights, m_tmp_flavour_weights);
    }
    inline void ResetFlavourSelectionWeights() {
      ResetTmpWeights(m_tmp_flavour_weights);
    }

    // Gluon splitting z selection (analytic density ratio, accepted z only)
    void GluonSplittingReweighting(const std::vector<double>& probs);
    inline void AcceptGluonSplittingWeights() {
      AcceptTmpWeights(m_gluon_weights, m_tmp_gluon_weights);
    }
    inline void ResetGluonSplittingWeights() {
      ResetTmpWeights(m_tmp_gluon_weights);
    }

    // Cluster splitting z selection (accept/reject reweighting per z trial)
    inline bool DoClusterSplittingReweighting(const int nsplit) const {
      return m_n_variations > 1 &&
             (m_reweight_max_nsplit < 0 || nsplit < m_reweight_max_nsplit);
    }
    void ClusterSplittingReweighting(const bool accepted,
                                     const std::vector<double>& probs);
    inline void AcceptClusterSplittingWeights() {
      AcceptTmpWeights(m_cluster_weights, m_tmp_cluster_weights);
    }
    inline void ResetClusterSplittingWeights() {
      ResetTmpWeights(m_tmp_cluster_weights);
    }

    // Transverse momentum selection (committed directly)
    void KTSelectionReweighting(const std::vector<double>& probs);

    // Hadron pair selection in soft cluster decays (committed directly)
    void SoftClusterReweighting(const std::vector<double>& weights,
                                const std::vector<double>& totweights);

    // Bernoulli prompt-decay decision, only stochastic for zeta > 0
    // (committed directly)
    void PromptDecayReweighting(const bool decayed,
                                const std::vector<double>& probs);

    void ApplyVariationWeights(ATOOLS::Blob * blob);
  };
}

#endif
