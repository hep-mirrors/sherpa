#ifndef AHADIC_Tools_Hadronisation_Reweighting_H
#define AHADIC_Tools_Hadronisation_Reweighting_H

#include "ATOOLS/Phys/Blob.H"
#include <vector>
#include <map>
#include <string>
// OUTPUT
#include "ATOOLS/Phys/Flavour.H"
#include <utility>
#include <fstream>

namespace ATOOLS { class Histogram; }

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

    // OUTPUT
    int  m_output_mode;
    bool m_reweighting_output;
    bool m_event_weights_applied;
    size_t m_total_events;
    std::vector<double> m_event_variation_weights;
    std::vector<double> m_event_flavour_weights;
    std::vector<double> m_event_gluon_weights;
    std::vector<double> m_event_cluster_weights;
    std::vector<double> m_event_soft_weights;
    std::vector<double> m_event_kt_weights;
    std::vector<double> m_event_promptdecay_weights;
    std::vector<double> m_sum_weights;
    std::vector<double> m_sum_weights_squared;
    std::vector<size_t> m_cutoff_count;
    void AccumulateEventStatistics();
    void PrintVariationStatistics();
    std::ofstream m_hadronisation_weight_file;
    typedef std::vector<std::pair<char, std::string> > Record_Vector;
    Record_Vector m_tmp_flavour_records;
    Record_Vector m_tmp_gluon_records;
    Record_Vector m_tmp_cluster_records;
    Record_Vector m_call_records;
    Record_Vector m_event_records;
    void CommitTmpRecords(Record_Vector& tmp_records);
    void DiscardTmpRecords(Record_Vector& tmp_records);
    std::string m_hist_path;
    std::map<std::string, ATOOLS::Histogram*> m_hist_nominal;
    std::map<std::string, std::vector<ATOOLS::Histogram*> > m_hist_reweighted;
    void BookHistograms();
    void FillHistograms();
    void WriteHistograms();
    void FillObs(const std::string& key, const double value);
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
      CommitTmpRecords(m_tmp_flavour_records); // OUTPUT
    }
    inline void ResetFlavourSelectionWeights() {
      ResetTmpWeights(m_tmp_flavour_weights);
      DiscardTmpRecords(m_tmp_flavour_records); // OUTPUT
    }

    // Gluon splitting z selection (analytic density ratio, accepted z only)
    void GluonSplittingReweighting(const std::vector<double>& probs);
    inline void AcceptGluonSplittingWeights() {
      AcceptTmpWeights(m_gluon_weights, m_tmp_gluon_weights);
      CommitTmpRecords(m_tmp_gluon_records); // OUTPUT
    }
    inline void ResetGluonSplittingWeights() {
      ResetTmpWeights(m_tmp_gluon_weights);
      DiscardTmpRecords(m_tmp_gluon_records); // OUTPUT
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
      CommitTmpRecords(m_tmp_cluster_records); // OUTPUT
    }
    inline void ResetClusterSplittingWeights() {
      ResetTmpWeights(m_tmp_cluster_weights);
      DiscardTmpRecords(m_tmp_cluster_records); // OUTPUT
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

    // OUTPUT
    void ResetStats();
    void WriteEventStatistics();
    void RecordFlavourPop(const ATOOLS::Flavour& flav, const double& Emax);
    void RecordGluonZ(const double& z);
    void RecordClusterZ(const int nsplit,
                        const int type1, const double& z1,
                        const int type2, const double& z2,
                        const double& mass);
    void RecordKT(const int site, const double& kt, const double& ktmax);
    void RecordSoftDecay(const double& mass, const ATOOLS::Flavour& had1,
                         const ATOOLS::Flavour& had2);
    void RecordTransition(const double& mass, const ATOOLS::Flavour& had);
    void RecordPrimaryCluster(const double& mass, const ATOOLS::Flavour& fl1,
                              const ATOOLS::Flavour& fl2, const bool direct);
  };
}

#endif
