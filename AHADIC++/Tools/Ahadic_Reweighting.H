#ifndef AHADIC_Tools_Ahadic_Reweighting_H
#define AHADIC_Tools_Ahadic_Reweighting_H

#include "ATOOLS/Phys/Blob.H"
#include <vector>
#include <map>
#include <string>

namespace AHADIC {

  class Ahadic_Reweighting {
  private:
    size_t m_n_variations;
    double m_max_reweight_factor;
    int    m_reweight_max_nsplit;
    std::map<std::string, std::vector<double>> m_variation_vectors;

    std::vector<double> m_variation_weights;
    // Per-channel event weights, accumulated over all splittings of an event:
    // flavour popping (Flavour_Selector), gluon splitting z (Gluon_Splitter),
    // cluster splitting z (Cluster_Splitter), soft cluster decays
    // (Soft_Cluster_Handler), and transverse momentum (KT_Selector).
    std::vector<double> m_flavour_weights;
    std::vector<double> m_gluon_weights;
    std::vector<double> m_cluster_weights;
    std::vector<double> m_soft_weights;
    std::vector<double> m_kt_weights;
    // Per-splitting weights, committed into the event weights when a
    // splitting succeeds and discarded when an attempt fails.
    std::vector<double> m_tmp_flavour_weights;
    std::vector<double> m_tmp_gluon_weights;
    std::vector<double> m_tmp_cluster_weights;
    std::vector<double> m_tmp_kt_weights;

    void AcceptTmpWeights(std::vector<double>& weights,
                       std::vector<double>& tmp_weights);
    void ResetTmpWeights(std::vector<double>& tmp_weights);
  public:
    Ahadic_Reweighting();
    ~Ahadic_Reweighting();

    void Initialize();
    void ResetEvent();

    inline bool   Active() const { return m_n_variations > 1; }
    inline size_t NumberOfVariations() const { return m_n_variations; }
    std::vector<double> GetVariationVector(const std::string& keyword) const;

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

    // Transverse momentum selection (accept/reject reweighting per kt trial)
    void KTSelectionReweighting(const bool accepted,
                                const std::vector<double>& probs);
    inline void AcceptKTSelectionWeights() {
      AcceptTmpWeights(m_kt_weights, m_tmp_kt_weights);
    }
    inline void ResetKTSelectionWeights() {
      ResetTmpWeights(m_tmp_kt_weights);
    }

    // Hadron pair selection in soft cluster decays (committed directly)
    void SoftClusterReweighting(const std::vector<double>& weights,
                                const std::vector<double>& totweights);

    void ApplyVariationWeights(ATOOLS::Blob * blob);
  };
}

#endif
